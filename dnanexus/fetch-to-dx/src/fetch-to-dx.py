#!/usr/bin/env python2.7
# fastqc 0.0.1
# Generated by dx-app-wizard.
#
# Parallelized execution pattern: Your app will generate multiple jobs
# to perform some computation in parallel, followed by a final
# "postprocess" stage that will perform any additional computations as
# necessary.
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import sys, os, subprocess, json, requests, shlex, urlparse, logging
from datetime import datetime
import dxpy

logger = logging.getLogger("Applet")

@dxpy.entry_point("postprocess")
def postprocess(reports):
    # Change the following to process whatever input this stage
    # receives.  You may also want to copy and paste the logic to download
    # and upload files here as well if this stage receives file input
    # and/or makes file output.

    #for output in reports:
    #   pass
    return reports

@dxpy.entry_point("noop")
def noop(enc_file_name, bucket_url, proj_id, dx_folder, file_acc, dx_file_name, skipvalidate=False):
    return {
        "report": None,
        "summary": None,
        "zip": None,
        "file": dx_file_name
    }

@dxpy.entry_point("process")
def process(enc_file_name, bucket_url, proj_id, dx_folder, file_acc, dx_file_name, skipvalidate=False):
    # Change the following to process whatever input this stage
    # receives.  You may also want to copy and paste the logic to download
    # and upload files here as well if this stage receives file input
    # and/or makes file output.

    print "* "+enc_file_name+" to "+dx_folder

    test = list( dxpy.find_data_objects(classname='file',
                           folder=dx_folder, project=proj_id, name_mode='exact',
                           name=dx_file_name, properties={ "accession": file_acc }, return_handler=False) )

    start = datetime.now()
    if not test or len(test) == 0:
        try:
            subprocess.check_call(shlex.split('aws s3 cp %s ./%s --quiet' %(bucket_url,dx_file_name)), stderr=subprocess.STDOUT)
            #subprocess.check_call(shlex.split('aws s3 cp %s ./%s' % (bucket_url,dx_file_name) ), stderr=subprocess.STDOUT)
        except:
            try:
                print "* s3 cp failed.  Reverting to 'wget'"
                web_url = "https://www.encodeproject.org/files/%s/@@download/%s" % (file_acc,enc_file_name)
                subprocess.check_call(shlex.split('wget %s -O %s --quiet' % (web_url,dx_file_name) ), stderr=subprocess.STDOUT)
            except:
                print "* ERROR: Upload failed"
                return {
                    "file": None,
                    "report": None,
                    "summary": None,
                    "zip": None
                }
        end = datetime.now()
        duration = end - start
        start = end
        print "* copied to dx local in %.2f seconds" % duration.seconds
 
        subprocess.check_call(shlex.split('ls -l %s' %(dx_file_name)))
        
        # Make sure folder exists before copying!
        project = dxpy.DXProject(proj_id)  ## should be default

        dx_file = dxpy.upload_local_file(dx_file_name, project=proj_id, folder=dx_folder, properties={ "accession": file_acc })
        end = datetime.now()
        duration = end - start
        print "* Uploaded to dx project in %.2f seconds" % duration.seconds

    else:
        dxpy.download_dxfile(test[0]['id'], dx_file_name)
        dx_file=dxpy.dxfile.DXFile(test[0]['id'])
        end = datetime.now()
        duration = end - start
        print "* Downloaded already existing file from in %.2f seconds" % duration.seconds

    if skipvalidate or not (dx_file_name.endswith(".fastq.gz") or dx_file_name.endswith(".fq.gz")):
        return {
            "file": dx_file,
            "report": None,
            "summary": None,
            "zip": None
        }

    subprocess.check_call(['mkdir', 'output'])
    print "* Run QC"
    fqc_command = "/usr/bin/FastQC/fastqc " + dx_file_name + " -o output"
    print "* " + fqc_command
    subprocess.check_output(shlex.split(fqc_command))
    subprocess.check_output(['ls','-l', 'output'])
    reads_basename = dx_file_name.rstrip('.gz').rstrip('.fq').rstrip('.fastq')
    subprocess.check_call(['unzip', "output/%s_fastqc.zip" % reads_basename])
    print "* Upload results"

    subprocess.check_call(['mv', "%s_fastqc/fastqc_data.txt" % reads_basename, "%s_data.txt" % reads_basename ])
    subprocess.check_call(['mv', "%s_fastqc/summary.txt" % reads_basename, "%s_summary.txt" % reads_basename ])

    report_dxfile = dxpy.upload_local_file("%s_data.txt" % reads_basename, folder=dx_folder, project=proj_id)
    summary_dxfile = dxpy.upload_local_file("%s_summary.txt" % reads_basename, folder=dx_folder, project=proj_id)
    zip_dxfile = dxpy.upload_local_file("output/%s_fastqc.zip" % reads_basename, folder=dx_folder, project=proj_id)
    print report_dxfile
    return {
        "file": dx_file,
        "report": report_dxfile,
        "summary": summary_dxfile,
        "zip": zip_dxfile
    }

@dxpy.entry_point("main")
def main(exp_acc, files_to_fetch=None, skipvalidate=True, key='www', debug=False):

    # Splits the work into parallel tasks: one for each file to fetch.

    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    proj_id = os.environ['DX_PROJECT_CONTEXT_ID']
    project = dxpy.DXProject(proj_id)  ## should be default

    logger.debug("* Project: " + proj_id)

    if files_to_fetch != None:
        logger.debug("* f2f_json: " + files_to_fetch)
        file_objs = json.loads(files_to_fetch.encode('ascii')) # Expect [ {},{},{},... ]
        logger.debug(file_objs)
            
    # f_obj = { "accession": ,"dx_folder": ,"dx_file_name": ,"enc_file_name": ,"bucket_url": }
                 
    subjobs = []
    if file_objs:
        for f_obj in file_objs:
            skipvalidate_this = skipvalidate
            dx_file_name = f_obj["dx_file_name"]
            if dx_file_name.endswith(".fastq.gz") or dx_file_name.endswith(".fq.gz"):
                skipvalidate_this = True

            logger.debug(f_obj["bucket_url"] + " " + f_obj["enc_file_name"])

            #process(f_obj["enc_file_name"], f_obj["bucket_url"], project.get_id(), f_obj["dx_folder"], f_obj["accession"], \
            #                                                                    f_obj["dx_file_name"], skipvalidate_this)
            subjob_input = {
                "enc_file_name": f_obj["enc_file_name"],
                "bucket_url": f_obj["bucket_url"],
                "proj_id": project.get_id(),
                "dx_folder": f_obj["dx_folder"],
                "file_acc": f_obj["accession"],
                "dx_file_name": f_obj["dx_file_name"],
                "skipvalidate": skipvalidate_this
            }
            subjobs.append(dxpy.new_dxjob(subjob_input, "process"))
            #subjobs.append(dxpy.new_dxjob(subjob_input, "noop"))

    # This does not wait for subjob completion as I thought.
    files_fetched = [subjob.get_output_ref("file") for subjob in subjobs]
    logger.debug("Attempting to fetch %d file(s)" % (len(files_fetched)))
    
    if skipvalidate:
        output = {
                    "fetched_count": len(files_fetched),
                    "files": files_fetched
        }
    else:
        output = {
                    "fetched_count": len(files_fetched),
                    "files": files_fetched,
                    "reports": [subjob.get_output_ref("report") for subjob in subjobs],
                    "summaries": [subjob.get_output_ref("summary") for subjob in subjobs],
                    "zips": [subjob.get_output_ref("zip") for subjob in subjobs],
        }


    return output

dxpy.run()
