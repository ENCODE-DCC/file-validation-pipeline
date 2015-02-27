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
                    required=False)

    ap.add_argument('-e', '--experiment',
                    help='ENCSR experiment accession',
                    required=False)

    ap.add_argument('--project',
                    help="Project to run analysis in (default: '" + PROJECT_DEFAULT + "')",
                    default=PROJECT_DEFAULT,
                    required=False)

    ap.add_argument('--assay',
                    help="Assay to loop over for all exps",
                    required=False)

    ap.add_argument('-t', '--all',
                    help='All replicates for assay type',
                    action='store_true',
                    required=False)

    ap.add_argument('--resultsLoc',
                    help="The location to to place results folders (default: '<project>:" + \
                                                                    RESULT_FOLDER_DEFAULT + "')",
                    default=RESULT_FOLDER_DEFAULT,
                    required=False)

    return ap

def get_fastqc(accession, project):
    summary_fn = accession+"_summary.txt"
    report_fn = accession+"_data.txt"

    summary_link = dxencode.find_file(summary_fn, project.get_id())
    report_link = dxencode.find_file(report_fn, project.get_id())

    metrics = {}
    try:
        with dxpy.open_dxfile(report_link) as rfd:
            total = re.compile('Total Sequences\s+(\d+)')
            for line in rfd:
                m = total.match(line)
                if m:
                    metrics.update({ 'Total Sequences': m.group(1) })
    except Exception, e:
        print "ERROR: Could not read FastQC summary: %s (%s) \n%s" % (summary_fn, summary_link, e)
        metrics.update({'Total Sequences': -999.999 })

    try:
        with dxpy.open_dxfile(summary_link) as sfd:
            fastqc = re.compile('(PASS|FAIL|WARN)\s+(.+)\s+ENCFF')
            for line in sfd:
                m = fastqc.match(line)
                if m:
                    metrics.update({ m.group(2):  m.group(1) })

    except Exception, e:
        print "ERROR: Could not read FastQC report: %s (%s) \n%s" % (report_fn, report_link, e)

    #print json.dumps(metrics)
    return metrics

def get_analysis_time(accession, repstr, project):

    result = list(dxpy.find_analyses(project=project, name='*'+accession+repstr+'*', name_mode='glob', describe=True, state='done'))
    if len(result) != 1:
        print "WARN: No single (%s) analysis found for %s%s" % (len(result), accession, repstr)
        return -999.99

    anl = result[0]['describe']
    start = anl['created']
    finish = [ t['setAt'] for t in anl['stateTransitions'] if t['newState'] == 'done'][0]

    #print start, finish, finish-start, (finish-start)/1000.0
    return (finish-start)/1000.0 # covert to secs.

def get_exp_time(accession, project, skip=False):

        expr = dxencode.encoded_get(SERVER+accession, AUTHID=AUTHID, AUTHPW=AUTHPW)
        try:
            expr.raise_for_status()
        except:
            print "ERROR: Could not find %s in db" % accession
            return

        exp = expr.json()

        if skip:
            return exp

        elapsed = {}
        for rep in exp['replicates']:
            repstr = '_rep'+str(rep['biological_replicate_number'])+'_'+str(rep['technical_replicate_number'])
            if not elapsed.get(repstr, None):
                elapsed[repstr] = { 'time': get_analysis_time(exp['accession'],repstr,project),
                                    'fastqs': []
                                  }

        fqs = [ f for f in exp['files'] if f['file_format'] == 'fastq']
        total_reads = {}
        paired = ""
        for fq in fqs:
            rep = fq['replicate']
            repstr = '_rep'+str(rep['biological_replicate_number'])+'_'+str(rep['technical_replicate_number'])
            total_reads[fq['accession']] = get_fastqc(fq['accession'], project)['Total Sequences']
            elapsed[repstr]['fastqs'].append(fq['accession'])
            if fq.get('paired_end', None):
                paired = 'Paired'

        for repstr in elapsed.keys():
            sizes = "\t".join([ "\t".join((fq, str(total_reads[fq]))) for fq in elapsed[repstr]['fastqs']])
            print "\t".join((exp['accession'],repstr, str(elapsed[repstr]['time']), sizes, paired))

        return exp

# woo hoo global
(AUTHID,AUTHPW,SERVER) = dxencode.processkey('default')

def main():
    argparser = get_args()
    args = argparser.parse_args()

    project = dxencode.get_project(args.project)

    if args.file:
        getr = dxencode.encoded_get(SERVER+args.file, AUTHID=AUTHID, AUTHPW=AUTHPW)
        try:
            getr.raise_for_status()
        except:
            print "Could not find %s in db" % args.file
            raise
        encff = getr.json()
        metrics = get_fastqc(encff['accession'], project)
        print json.dumps(metrics)
    elif args.experiment:
        get_exp_time(args.experiment, project)
    elif args.all:
        assay = args.assay or "OBI:0001271"
        query = '/search/?type=experiment&assay_term_id=%s&award.rfa=ENCODE3&limit=all&frame=embedded&replicates.library.biosample.donor.organism.name=mouse&files.file_format=fastq' % assay

        res = dxencode.encoded_get(SERVER+query, AUTHID=AUTHID, AUTHPW=AUTHPW)
        exps = res.json()['@graph']

        for exp in exps:
            acc = exp['accession']
            if len(exp['replicates']) > 0:
                if exp['replicates'][0]['library'].get('size_range', "") != '>200':
                    print "Skipping %s with wrong library size (%s)" % (acc, exp['replicates'][0]['library'].get('size_range', ""))
                    #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                    continue
                if exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity_units', "") == "cells":
                    ncells = float(exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity', 0.0))
                    if ncells < 20:
                        print "Skipping %s as single-cell (%s %s)" % (acc, exp['replicates'][0]['library'].get('nucleic_acid_starting_quantity_units', ""), ncells)
                        #print json.dumps(exp['replicates'][0]['library'], sort_keys=True, indent=4, separators=(',',': '))
                        continue
            get_exp_time(acc, project)


if __name__ == '__main__':
    main()
