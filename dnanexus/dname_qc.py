#!/usr/bin/env python

import dxpy
import argparse
import re
import json
from dxencode import dxencode as dxencode
import get_fastqc

PROJECT_DEFAULT = 'dna-me-pipeline'
RESULT_FOLDER_DEFAULT = '/runs'

APPROX_HAPLOID_GENOME_SIZE = 3e09
## WARNING HACK

labels = [
        'Sequences analysed in total',
        'Mapping efficiency',
        'C methylated in CpG context',
        'C methylated in CHG context',
        'C methylated in CHH context'
]


def parse_map_report(folder, project):

    mapreport = "/*_bismark_map_report.txt"
    report_link = dxencode.find_file(folder+mapreport, project.get_id(), recurse=False)

    metrics = {}
    res = {}
    for lab in labels:
        res[lab] = re.compile("(%s):\s+(.+)" % lab)

    try:
        with dxpy.open_dxfile(report_link) as rfd:
            for line in rfd:
                m = False
                for metric in res.values():
                    m = metric.match(line)
                    if m:
                        metrics.update({ m.group(1): m.group(2).strip() })
                        continue

    except Exception, e:
        print "ERROR: Could not read Bismark mapping report in %s (%s) \n%s" % (folder, report_link, e)

    return metrics

def get_bismark_stats(experiment, project):
    mbias = "*_bismark.M-bias.txt"

    metrics = {}

    for rep in experiment['replicates']:
        acc = experiment['accession']
        metrics[acc] = {}

        rep_str = "%s_%s" % (rep['biological_replicate_number'], rep['technical_replicate_number'])
        folder = RESULT_FOLDER_DEFAULT+'/'+acc+'/'+rep_str

        mapped = parse_map_report(folder, project)
        lambda_mapped = parse_map_report(folder+'/lambda', project)

        nreads = float(mapped.get('Sequences analysed in total', -999))
        mapeff = float(mapped.get('Mapping efficiency', '-99999%').strip('%'))/100.0
        coverage = nreads * mapeff * float(rep['read_length']) / APPROX_HAPLOID_GENOME_SIZE
        print '\t'.join([acc,rep_str]+[ mapped.get(v,'-999.9') for v in labels ] + [ lambda_mapped.get(v,'-999.9') for v in labels ]+["%2.2f" % coverage])

(AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')

def process_exp(acc, project, skipfq):

    fqc_metrics = {}

    expr = get_fastqc.get_exp_time(acc, project, skip=skipfq)
    if not skipfq:
        fqc_metrics[expr['accession']] = {}
        for fq in [ f for f in expr['files'] if f['file_format'] == 'fastq' ]:
            fqc_metrics[expr['accession']][fq['accession']] = get_fastqc.get_fastqc(fq['accession'], project)
    #print json.dumps(fqc_metrics, indent=4)
    get_bismark_stats(expr, project)

def main():
    argparser = get_fastqc.get_args()

    argparser.add_argument('--skipfq',
                            help='Skip parsing fastqc',
                            action='store_true',
                            required=False)

    args = argparser.parse_args()

    project = dxencode.get_project(args.project)
    print "\t".join(['Experiment','Replicate']+labels+["lambda "+l for l in labels]+['Estimated Coverage'])
    if args.experiment:
        process_exp(args.experiment, project, skipfq=args.skipfq)
    elif args.all:
        assay = args.assay or "OBI:0001863"
        query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&files.file_format=fastq' % assay

        res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
        exps = res.json()['@graph']

        for exp in exps:
            acc = exp['accession']
            if len(exp['replicates']) > 0:
               process_exp(acc, project, skipfq=args.skipfq)


if __name__ == '__main__':
    main()
