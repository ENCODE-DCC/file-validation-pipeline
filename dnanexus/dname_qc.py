#!/usr/bin/env python

import dxpy
import argparse
import re
import json
from dxencode import dxencode as dxencode
import get_fastqc

PROJECT_DEFAULT = 'dna-me-pipeline'
RESULT_FOLDER_DEFAULT = '/runs'


def parse_map_report(folder, project):

    mapreport = "*_bismark_map_report.txt"
    report_link = dxencode.find_file(mapreport, project.get_id(), folder=folder, recurse=False)

    metrics = {}
    labels = [
        'Sequences analysed in total'
        'Mapping efficiency',
        'C methylated in CpG context',
        'C methylated in CHG context',
        'C methylated in CHH context'
    ]
    res = {}
    for lab in labels:
        res[lab] = re.compile("(%s):\s+(.+)\s")

    try:
        with dxpy.open_dxfile(report_link) as rfd:
            for line in rfd:
                for metric in res.keys():
                    m = metric.match(line)
                if m:
                    metrics.update({ m.group(1): m.group(2) })
        rfd.close()

    except Exception, e:
        print "ERROR: Could not read Bismark mapping report in %s (%s) \n%s" % (folder, report_link, e)


def get_bismark_stats(experiment, project):
    mbias = "*_bismark.M-bias.txt"

    metrics = {}

    for rep in experiment['replicates']:
        acc = experiment['accession']
        metrics[acc] = {}

        rep_str = rep['biological_replicate_number']+'_'+rep['technical_replicate_number']
        folder = RESULT_FOLDER_DEFAULT+'/'+acc+'/'+rep_str

        metrics[acc][rep_str].extend(get_qc_from_reports(folder))
    return metrics

(AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')

def main():
    args = get_args()

    project = dxencode.get_project(args.project)
    fqc_metrics = {}

    if args.experiment:
        expr = get_fastqc.get_exp_time(args.experiment, project)
        for fq in [ f for f in expr['files'] if f['file_format'] == 'fastq' ]:
            fqc_metrics[(expr['accession'],fq['accession'])] = get_fastqc.get_fastqc(fq['accession'], project)
        print json.dumps(fqc_metrics, indent=4)
        get_bismark_stats(expr, project)
    elif args.all:
        assay = args.assay or "OBI:0001863"
        query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&files.file_format=fastq' % assay

        res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
        exps = res.json()['@graph']

        for exp in exps:
            acc = exp['accession']
            if len(exp['replicates']) > 0:
                get_fastqc.get_exp_time(acc, project)


if __name__ == '__main__':
    main()
