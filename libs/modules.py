import os
import logging
import subprocess


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


def run_cutadapt(local_path, s_out, sample_id, s_input):
    logger = get_logger()
    logger.info('Start run cutadapt: {}'.format(sample_id))
    execute = os.path.join(local_path, 'bin/TrimGalore-0.6.5/trim_galore')
    cutadapt_path = os.path.join(local_path, 'bin/cutadapt/cutadapt_program')
    outdir = os.path.join(s_out, 'clean_data')
    fq1 = os.path.join(s_input, '{}_1.fq.gz'.format(sample_id))
    fq2 = os.path.join(s_input, '{}_2.fq.gz'.format(sample_id))
    if os.path.exists(fq1) and os.path.exists(fq2):
        pass
    else:
        exit(1)
    os.system('mkdir -p {}'.format(outdir))

    cmd = '{} ' \
          '--paired ' \
          '--path_to_cutadapt {} ' \
            '--stringency 3 ' \
            '-j 4 ' \
            '--gzip ' \
            '--trim-n ' \
            '--length 75 ' \
            '-o {} ' \
            '--no_report_file ' \
            '-q 30 ' \
            '--basename {} ' \
            '{} ' \
            '{} >> {}/cutadapt.log 2>&1 && touch {}/{}.SUCCESS'.format(execute, cutadapt_path, outdir, sample_id, fq1, fq2, outdir, outdir, sample_id)
    # print(cmd)
    if os.path.exists('{}/{}.SUCCESS'.format(outdir, sample_id)):
        logger.info('{} cutadapt finished.Skipped.'.format(sample_id))
        return True
    else:
        pass
    cutadatp = subprocess.Popen(cmd, shell=True)
    cutadatp.wait()
    if os.path.exists('{}/{}.SUCCESS'.format(outdir, sample_id)):
        return True
    else:
        return False
    # logger.info('Cutadapt done')


def run_hisat2(local_path, s_out, sample_id, index):
    logger = get_logger()
    logger.info('Start run hisat2: {}'.format(sample_id))
    hisat2_execute = os.path.join(local_path, 'bin/hisat2-2.1.0/hisat2')
    samtools_execute = os.path.join(local_path, 'bin/samtools-1.10/samtools')
    outdir = os.path.join(s_out, 'alignment_bams')
    fq1 = os.path.join(s_out, 'clean_data', '{}_val_1.fq.gz').format(sample_id)
    fq2 = os.path.join(s_out, 'clean_data', '{}_val_2.fq.gz').format(sample_id)
    if os.path.exists(fq1) and os.path.exists(fq2):
        pass
    else:
        exit(1)
    os.system('mkdir -p {}'.format(outdir))
    cmd = '{} ' \
          '-x {} ' \
          '-p 4 ' \
          '-1 {} ' \
          '-2 {} | ' \
          '{} view ' \
          '-@ 4 ' \
          '-Sb | ' \
          '{} sort ' \
          '-@ 4 ' \
          '-m 2G ' \
          '-o {}/{}_sorted.bam && ' \
          '{} index ' \
          '-@ 4 ' \
          '{}/{}_sorted.bam ' \
          '{}/{}_sorted.bam.bai > {}/align.log 2>&1 && touch {}/{}.SUCCESS'.format(hisat2_execute, index, fq1, fq2, samtools_execute, samtools_execute, outdir,
                                        sample_id, samtools_execute, outdir, sample_id, outdir, sample_id, outdir, outdir, sample_id)
    # print(cmd)
    if os.path.exists('{}/{}.SUCCESS'.format(outdir, sample_id)):
        logger.info('{} alignment finished. Skipped.'.format(sample_id))
        return True
    else:
        pass
    alignment = subprocess.Popen(cmd, shell=True)
    alignment.wait()
    if os.path.exists('{}/{}.SUCCESS'.format(outdir, sample_id)):
        return True
    else:
        return False
    # logger.info('Hisat2 done')


def run_featurecount(local_path, s_out, species, logger):
    logger.info('Start running featureCount')
    execute = os.path.join(local_path, 'bin/subread-2.0.1-Linux-x86_64/bin/featureCounts')
    outdir = os.path.join(s_out, 'counts')
    bam_dir = os.path.join(s_out, 'alignment_bams')
    if os.path.exists(bam_dir):
        pass
    else:
        exit(1)
    gtf = species[list(species.keys())[0]]['gtf']
    os.system('mkdir -p {}'.format(outdir))
    cmd = '{} ' \
          '-T 8 ' \
          '-p ' \
          '-g gene_id ' \
          '-a {} ' \
          '-o {}/feature_count.txt ' \
          '{}/*.bam && touch {}/feature_count.SUCCESS'.format(execute, gtf, outdir, bam_dir, outdir)
    # print(cmd)
    if os.path.exists('{}/feature_count.SUCCESS'.format(outdir)):
        logger.info('Call counts finished. Skipped.')
        return True
    featurecount = subprocess.Popen(cmd, shell=True)
    featurecount.wait()
    if os.path.exists('{}/feature_count.SUCCESS'.format(outdir)):
        return True
    else:
        return False
    # logger.info('FeatureCount done')


def run_overall_plot(local_path, s_out, s_sample_list):
    logger = get_logger()
    logger.info('Start run overall_plot')
    s_out = os.path.abspath(s_out)
    counts = '{}/counts/tmp_count.txt'.format(s_out)
    output = '{}/DE_analysis'.format(s_out)
    output = os.path.abspath(output)
    try:
        os.mkdir(output)
    except:
        pass
    s_sample_list = os.path.abspath(s_sample_list)
    if os.path.exists(counts):
        pass
    else:
        logger.error('{} cannot be found. Exit with code 1'.format(counts))
        return False
    cmd = 'Rscript ' \
          '{}/libs/overall_plot.R ' \
          '{} ' \
          '{} ' \
          '{} > ' \
          '{}/overall_plot.log 2>&1 ' \
          '&& touch {}/overall_plot.SUCCESS'.format(local_path, counts, s_sample_list, output, output, output)
    # print(cmd)
    overall_plot = subprocess.Popen(cmd, shell=True)
    overall_plot.wait()
    if os.path.exists('{}/overall_plot.SUCCESS'.format(s_out)):
        logger.info('Overall plot already finished. Skipped')
        return True
    else:
        exit(1)
        # overall_plot = subprocess.Popen(cmd, shell=True)
        # overall_plot.wait()
        # if os.path.exists('{}/overall_plot.SUCCESS'.format(s_out)):
        #     return True
        # else:
        #     return False


def run_deseq2(local_path, s_out, s_sample_list, species, s_resoud, s_control, treat, logger, pvalue, padjust, lfc, plot_type):
    logger.info('Start run deseq2 and enrichment: {}vs{}'.format(treat, s_control))
    s_out = os.path.abspath(s_out)
    s_sample_list = os.path.abspath(s_sample_list)
    counts = '{}/counts/tmp_count.txt'.format(s_out)
    if os.path.exists(counts):
        pass
    else:
        logger.error('{} cannot be found. Exit with code 1'.format(counts))
        exit(1)
    s_resoud = os.path.abspath(s_resoud)
    os.system('mkdir -p {}/Enrichment'.format(s_resoud))
    cmd = 'Rscript ' \
          '{}/libs/DESeq2.R ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} ' \
          '{} > {}/DESeq2.log 2>&1 && ' \
          'touch {}/SUCCESS'.format(local_path, counts, s_sample_list, list(species.keys())[0], s_resoud, s_control, treat, pvalue, padjust, lfc , plot_type, s_resoud, s_resoud)
    # print(cmd)
    if os.path.exists('{}/SUCCESS'.format(s_resoud)):
        logger.info('{} compare finished. Skipped.'.format(treat + '_vs_' + s_control))
        return True
    else:
        pass
    DESeq2 = subprocess.Popen(cmd, shell=True)
    DESeq2.wait()
    try:
        os.remove('Rplots.pdf')
    except:
        return True
    if os.path.exists('{}/SUCCESS'.format(s_resoud)):
        return True
    else:
        return False
    # logger.info('DESeq2 and enrichment done')
