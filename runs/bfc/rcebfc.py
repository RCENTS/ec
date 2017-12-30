from SCons.Script import Command
from SCons.Script import Environment


def sra_prefetch_cmd(data_dir, acc):
    pfetch_cmd = 'prefetch {}'
    pfetch_run = pfetch_cmd.format(acc)
    pfetch_tgt = "{}/sra/{}.sra".format(data_dir, acc)
    return Command(pfetch_tgt, [], pfetch_run)


def sra_validate_cmd(data_dir, acc, pfetch = []):
    validate_cmd = 'vdb-validate {}'
    validate_run = validate_cmd.format(acc)
    valid_tgt = "{}/sra/{}.vld".format(data_dir, acc)
    return Command(valid_tgt, pfetch, validate_run + " 2> $TARGET")

def sra_fastq_dump_cmd(data_dir, acc, paired=True, valid = []):
    fdump_cmd = 'fastq-dump -O {}/fastq/ -I --split-files {}'
    fdump_run = fdump_cmd.format(data_dir, acc)
    fdump_tgts = [ "{}/fastq/{}_1.fastq".format(data_dir, acc),
                   "{}/fastq/{}_2.fastq".format(data_dir, acc)]
    return Command(fdump_tgts, valid, fdump_run)


def sra_download_cmd(data_dir, acc, paired=True):
    return sra_fastq_dump_cmd(data_dir, acc, paired,
                              sra_validate_cmd(data_dir, acc,
                                               sra_prefetch_cmd(data_dir, acc)))
def ec_bfc_cmd(data_dir, acc, params):
    fdump_out = [ "{}/fastq/{}_1.fastq".format(data_dir, acc),
                   "{}/fastq/{}_2.fastq".format(data_dir, acc)]
    bfc_cmd = 'bfc {} {}'
    bfc_input = "{}/fastq/{}.fastq".format(data_dir, acc)
    bfc_output = "{}/fastq/{}_bfc.fastq".format(data_dir, acc)
    return Command(bfc_output,
                   Command(bfc_input, fdump_out, "cat $SOURCES > $TARGET"),
                   bfc_cmd.format(params, bfc_input) + " > $TARGET")


rdwn = sra_download_cmd('/data', 'SRR065390')
rbfc = ec_bfc_cmd('/data', 'SRR065390' , '-t8 -s 100m -k 23')
