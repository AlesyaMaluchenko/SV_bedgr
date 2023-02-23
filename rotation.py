import io
import os
import pandas as pd
from subprocess import call, check_output
import argparse
import multiprocessing as mp


#function takes vcf file and returns three BedGraph files with different SV (dup, del, inv) (filtered by quality)
def splitter(path_to_file, output, qual):
    file_name = os.path.basename(path_to_file)
    with open(path_to_file, 'r') as file:
        res = []
        for line in file:
            if line[0] == '#' and line[1] == '#':
                pass
            else:
                res.append(line)
        pd_file = pd.read_csv(io.StringIO(''.join(res)), sep='\t')
        end = []
        for i in range(len(pd_file)):
            end.append(pd_file['INFO'][i].split(';')[2][4:])
        pd_file['END'] = end
        pd_file['NUMB'] = 1
        for pref in ['DUP', 'DEL', 'INV']:
            spl = pd_file.loc[pd_file['ALT'] == '<{}>'.format(pref)]
            if qual:
                spl = spl.loc[spl['QUAL'] >= qual]
            name = file_name[0: -3] + '{}.BedGraph'.format(pref)
            bed = spl.to_csv(os.path.join(output, name), columns=['#CHROM', 'POS', 'END', 'NUMB'], sep='\t', header=False, index=False)
            with open(os.path.join(output, name), 'r') as file:
                data = file.read()
                with open(os.path.join(output, name), 'w') as new_file:
                    new_file.write("track type=bedGraph name=\"{pref} from {file}\" description=\"{pref} with quality > {quality}\" color=100,100,100\n".format(pref=pref, file=file_name, quality=qual) + data)


#calculation of frequency of SV in each interval
def frec_calc(x):
    return round((len(x) - x.tolist().count(0.0)) / len(x), 4)

def freq_calc_main(out, pref, filter_thresh, bedtools): 
    table = pd.read_table(os.path.join(out, '{}_combined.BedGraph'.format(pref)), header=None)
    selection = table.iloc[:, 0:3]
    selection['freq'] = table.iloc[0:, 3:].apply(frec_calc, axis=1)
    selection.to_csv(os.path.join(out, '{}_FREQ.BedGraph'.format(pref)), sep="\t", index=False, header=False)
    selection.query("freq > {}".format(filter_thresh)).to_csv(
        os.path.join(out, '{}_FREQ_filtered.BedGraph'.format(pref)), sep="\t", index=False, header=False)
    with open(os.path.join(out, '{}_FREQ_filtered.BedGraph'.format(pref)), 'r') as file:
        data = file.read()
        with open(os.path.join(out, '{}_FREQ_filtered.BedGraph'.format(pref)), 'w') as new_file:
            new_file.write("track type=bedGraph name=\"{pref}\" description=\"{pref} with frequency > {frequency}\" color=100,100,100\n".format(pref=pref, frequency=filter_thresh) + data)


def call_bash(command_bash):
    call(command_bash, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('vcfs_dir', nargs='?',
                        help='path to folder containing vcf files')
    parser.add_argument('bedgr_dir', nargs='?', help='path to folder containing BedGraph files after splitting')
    parser.add_argument('bedtools', nargs='?', help='path to bedtools exec', default='/usr/bin/bedtools')
    parser.add_argument('output', nargs='?', help='path to output folder')
    parser.add_argument('--threads', nargs='?', default=50, help='number of parallel threads for splitting files')
    parser.add_argument('--quality', nargs='?', default=0, help='Intervals with quality above X will be used')
    parser.add_argument('--filter_thresh', nargs='?',
                        default=0.01, help='Intervals with freq above X will be used for filtering in the pipeline')

    args = parser.parse_args()

    if not os.path.isdir(args.bedgr_dir):
        call('mkdir {}'.format(args.bedgr_dir), shell=True)

    if not os.path.isdir(args.output):
        call('mkdir {}'.format(args.output), shell=True)

    print('Start splitting VCFs')
    f_list = os.listdir(args.vcfs_dir)
    pool = mp.Pool(int(args.threads))
    for vcf in f_list:
        if vcf[-3:] == 'vcf':
            pool.apply_async(splitter, args=(os.path.join(args.vcfs_dir, vcf), args.bedgr_dir, args.quality))
    pool.close()
    pool.join()
    print('Splitting done')
    

    pool = mp.Pool(4)
    print('Union of separate tracks into single matrix')
    for pref in ['DUP', 'DEL', 'INV']:
        commandline = "{bedtools} unionbedg -i {bedgr_dir}/*{pref}.BedGraph > {out}".format(
            bedgr_dir=args.bedgr_dir,
            out=os.path.join(args.output, '{}_combined.BedGraph'.format(pref)),
            pref=pref,
            bedtools=args.bedtools)
        pool.apply_async(call_bash, args=(commandline,))
    pool.close()
    pool.join()
    print('Union done')

    pool = mp.Pool(4)
    print('Frequency calculation')
    for pref in ['DUP', 'DEL', 'INV']:
        pool.apply_async(freq_calc_main, args=(args.output, pref, args.filter_thresh, args.bedtools))
    pool.close()
    pool.join()
    print('Frequency calc done')

    call('{bedtools} unionbedg -i {DUP} {DEL} {INV} > {out}'.format(
        DUP=os.path.join(args.output, 'DUP_FREQ_filtered.BedGraph'),
        DEL=os.path.join(args.output, 'DEL_FREQ_filtered.BedGraph'),
        INV=os.path.join(args.output, 'INV_FREQ_filtered.BedGraph'),
        out=os.path.join(args.output, 'SV_combined.BedGraph'),
        bedtools=args.bedtools), shell=True)


    table = pd.read_table(os.path.join(args.output, 'SV_combined.BedGraph'.format(pref)), header=None)
    selection = table.iloc[:, 0:3]
    selection['max'] = table.iloc[0:, 3:].apply(max, axis=1)
    selection.to_csv(os.path.join(args.output, 'SV_combined.BedGraph'), sep="\t", index=False, header=False)



if __name__ == '__main__':
    main()



