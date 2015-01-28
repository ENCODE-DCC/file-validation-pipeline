#!/usr/bin/env python

import dxpy
import argparse
import re
import json
from dxencode import dxencode as dxencode

PROJECT_DEFAULT = 'long-rna-seq-pipeline'
RESULT_FOLDER_DEFAULT = '/runs'

def get_args():
    '''Parse the input arguments.'''
    ### LRNA specific
    ap = argparse.ArgumentParser(description="Gets FASTC results from dnanexus, creates QC object posts to encoded")
    ### LRNA specific

    ap.add_argument('-f', '--file',
                    help='ENCFF file accession',
                    required=True)

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    default=PROJECT_DEFAULT,
                    required=False)


    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    return ap.parse_args()

def main():
    args = get_args()

    (AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')
    print args
    getr = dxencode.encoded_get(SERVER+args.file, AUTHID=AUTHID, AUTHPW=AUTHPW)
    try:
        getr.raise_for_status()
    except:
        print "Could not find %s in db" % args.file
        raise
    encff = getr.json()
    summary_fn = encff['accession']+"_summary.txt"
    report_fn = encff['accession']+"_data.txt"

    project = dxencode.get_project(args.project)
    summary_link = dxencode.find_file(summary_fn, project.get_id())
    report_link = dxencode.find_file(report_fn, project.get_id())

    metrics = []
    with dxpy.open_dxfile(report_link) as rfd:
        total = re.compile('Total Sequences (\d+)')
        for line in rfd:
            m = total.match(line)
            if m:
                metrics.append({ 'metric': 'in total',
                                 'value':  m.group(0) })
    rfd.close()

    with dxpy.open_dxfile(summary_link) as sfd:
        fastqc = re.compile('(PASS|FAIL|WARN)\s+(.+)\s+ENCFF')
        for line in sfd:
            m = fastqc.match(line)
            if m:
                metrics.append({ 'metric': m.group(2),
                                 'value':  m.group(1) })
    sfd.close()

    print json.dumps(metrics)


if __name__ == '__main__':
    main()
