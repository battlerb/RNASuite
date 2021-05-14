import argparse
import logging
import sys
import multiprocessing
from libs import modules


def init_par():
    args = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                   description='''\
Author: Tianyu_Cui
Time: 21_11_2020
E-mail: 770210480@gmail.com
Version: v0.1(beta)
    ''',
                                   prog='RNA-Seq Pipeline (Basic Version)'
                                   )
    args.add_argument('-l', '--sample-list', help='输入样本分组信息', required=True, type=str)
    args.add_argument('-p', '--threads', help='程序并行数量 [默认为1, 参考：CPU核心数小于等于4、内存16G，设置并行数量为1；CPU核心数8、内存32G，设置并行数量为2]', default=1, type=int)
    args.add_argument('-o', '--output', help='输出结果目录 [默认为当然目录下的output文件夹]', type=str, default='output')
    args.add_argument('-i', '--input-dir', help='原始目录所在文件夹', required=True, type=str)
    args.add_argument('-s', '--species', help='物种类型', required=True, type=str)
    args.add_argument('-c', '--compare', help='基因差异表达分析比较文件', required=True, type=str)
    args.add_argument('-j', '--jump', help='跳过某几个步骤, {qc | align | count | de | overall_plot} [[0]]', type=str, default='0')
    args.add_argument('--pvalue', help='差异表达分析筛选的P-Value', default=0.05, type=float)
    args.add_argument('--padjust', help='差异表达分析筛选的P-Adjust', default=0.05, type=float)
    args.add_argument('--lfc', help='差异表达分析筛选的log2FoldChang绝对值', default=1.5, type=float)
    args.add_argument('--plot-type', help='输出图片的格式,默认为pdf和png都输出', default="both", type=str)
    args = args.parse_args()
    return args


def get_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    fmt = "%(asctime)-15s %(levelname)s %(filename)s %(lineno)d %(process)d %(message)s"
    datefmt = "%a %d %b %Y %H:%M:%S"
    formatter = logging.Formatter(fmt, datefmt)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def parse_sample_list(s_sample_list, logger):
    logger.info('Parsing sample list.')
    try:
        with open(s_sample_list, 'r') as sample_file:
            lines = sample_file.readlines()
            sample_ids = []
            sample_groups = []
            for i in range(1, len(lines)):
                sample_id = lines[i].strip().split(',')[0]
                sample_group = lines[i].strip().split(',')[-1]
                sample_ids.append(sample_id)
                sample_groups.append(sample_group)
            logger.info('Parsing sample list success')
            return sample_ids, sample_groups
    except:
        logger.critical('Cannot parse sample list file. Exit with code 1')
        exit(1)


def get_species(local_path, species, logger):
    species_info = {}
    species_data = {'hsa':
                        {'genome_index': '{}/db/Genome/GRCh38/Homo_sapiens.GRCh38'.format(local_path),
                         'gtf': '{}/db/Genome/GRCh38/Homo_sapiens.GRCh38.101.chr.gtf'.format(local_path)},
                    'mmu':
                        {'genome_index': '{}/db/Genome/GRCm38/Mus_musculus.GRCm38'.format(local_path),
                         'gtf': '{}/db/Genome/GRCm38/Mus_musculus.GRCm38.101.chr.gtf'.format(local_path)},
                    'ssc':
                        {'genome_index': '{}/db/Genome/Sscrofa11.1/Sus_scrofa.Sscrofa11.1'.format(local_path),
                         'gtf': '{}/db/Genome/GRCm38/Sus_scrofa.Sscrofa11.1.101.chr.gtf'.format(local_path)}}
    try:
        species_info[species] = species_data[species]
        logging.info('Your selected species is {}'.format(species))
    except KeyError:
        logger.warning('Only hsa and mmu are supported yet. Please wait for further update.')
        exit(1)
    return species_info


def parse_compare(s_compare, logger):
    logger.info('Start to parse compare file.')
    with open(s_compare, 'r') as compare_file:
        compare_dict = {}
        lines = compare_file.readlines()
        for i in range(1, len(lines)):
            compare_dict[i] = {'treat': lines[i].strip().split(',')[0], 'control': lines[i].strip().split(',')[1]}
    logger.info('Parse compare file done. You have {} group waiting to analyse'.format(len(compare_dict)))
    return compare_dict


def edit_count_file(s_out):
    original_count = '{}/counts/feature_count.txt'.format(s_out)
    with open(original_count, 'r') as original_count_file:
        lines = original_count_file.readlines()
    old_header = lines[1].strip().split('\t')
    new_header = []
    for i in range(len(old_header)):
        if i in range(6):
            new_header.append(old_header[i])
        else:
            new = old_header[i].split('/')[-1].replace('_sorted.bam', '')
            new_header.append(new)
    lines[1] = '\t'.join(new_header) + '\n'
    with open('{}/counts/tmp_count.txt'.format(s_out), 'w+') as new_count:
        for i in lines:
            new_count.write(i)


def run_analysis_1(s_out, s_sample_list, local_path, s_input, sample_ids, s_threads, species, groups, s_compare, logger, jump, pvalue, padjust, lfc, plot_type):
    jump = jump.split(',')
    pool1 = multiprocessing.Pool(processes=s_threads)
    for sample_id in sample_ids:
        if 'qc' not in jump or jump == ['0']:
            try:
                pool1.apply_async(modules.run_cutadapt, (local_path, s_out, sample_id, s_input))
            except:
                logger.error('Cannot running cutadapt')
                pass
        else:
            logger.info('Detect jump code "qc". Skip this step')
    pool1.close()
    pool1.join()

    pool1 = multiprocessing.Pool(processes=s_threads)
    for sample_id in sample_ids:
        if 'align' not in jump or jump == ['0']:
            try:
                index = species[list(species.keys())[0]]['genome_index']
                pool1.apply_async(modules.run_hisat2, (local_path, s_out, sample_id, index))
            except:
                logger.error('Cannot run hisat2')
                pass
        else:
            logger.info('Detect jump code "align". Skip this step')
    pool1.close()
    pool1.join()

    try:
        if 'count' not in jump or jump == ['0']:
            modules.run_featurecount(local_path, s_out, species, logger)
        else:
            logger.info('Detect jump code "count". Skip this step')
    except:
        logger.error('Cannot run featureCount')
        pass
    edit_count_file(s_out)
    try:
        if 'overall_plot' not in jump or jump == ['0']:
            modules.run_overall_plot(local_path, s_out, s_sample_list)
        else:
            logger.info('Detect jump code "overall_plot". Skip this step')
    except:
        logger.error('Cannot run overall_plot')
        pass
    pool2 = multiprocessing.Pool(processes=s_threads)
    if 'de' not in jump or jump == ['0']:
        for group in s_compare.keys():
            treat = s_compare[group]['treat']
            control = s_compare[group]['control']
            try:
                resoud = '{}/DE_analysis/{}_VS_{}'.format(s_out, treat, control)
                logger.info('Start run deseq2 and enrichment: {}vs{}'.format(treat, control))
                pool2.apply_async(modules.run_deseq2, (local_path, s_out, s_sample_list, species, resoud, control, treat, logger, pvalue, padjust, lfc, plot_type))
            except:
                logger.error('Cannot run DESeq2 and enrichment for {} group'.format(s_compare[group]))
    else:
        logger.info('Detect jump code "de". Skip this step')
    pool2.close()
    pool2.join()


def main(local_path):
    args = init_par()
    logger = get_logger()
    sys.path.append('{}/bin/coloredlogs_lib'.format(local_path))
    import coloredlogs
    coloredlogs.install(level='DEBUG')
    s_compare = parse_compare(args.compare, logger)
    sample_ids, sample_groups = parse_sample_list(args.sample_list, logger)
    species_info = get_species(local_path, args.species, logger)
    run_analysis_1(args.output, args.sample_list, local_path, args.input_dir,
                 sample_ids, args.threads, species_info, sample_groups, s_compare, logger, args.jump, args.pvalue, args.padjust, args.lfc, args.plot_type)
