#!/usr/bin/env python
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
import dxpy

KEYFILE = 'keypairs.json'
DEFAULT_SERVER = 'https://www.encodeproject.org'
S3_SERVER='s3://encode-files/'

logger = logging.getLogger("Applet")

def processkey(key):

    if key:
        keysf = open(KEYFILE,'r')
        keys_json_string = keysf.read()
        keysf.close()
        keys = json.loads(keys_json_string)
        key_dict = keys[key]
    else:
        key_dict = {}
    AUTHID = key_dict.get('key')
    AUTHPW = key_dict.get('secret')
    if key:
        SERVER = key_dict.get('server')
    else:
        SERVER = 'https://www.encodeproject.org/'

    if not SERVER.endswith("/"):
        SERVER += "/"

    return (AUTHID,AUTHPW,SERVER)

def encoded_get(url, AUTHID=None, AUTHPW=None):
    HEADERS = {'content-type': 'application/json'}
    if AUTHID and AUTHPW:
        response = requests.get(url, auth=(AUTHID,AUTHPW), headers=HEADERS)
    else:
        response = requests.get(url, headers=HEADERS)
    return response

def find_or_create_folder(project, sub_folder, root_folder='/'):
    folder = root_folder+sub_folder
    logger.debug("Creating %s (%s)" % (folder, root_folder))
    if folder in project.list_folder(root_folder)['folders']:
        return folder
    else:
        return project.new_folder(folder)

def get_bucket(SERVER, AUTHID, AUTHPW, f_obj):

    #make the URL that will get redirected - get it from the file object's href property
    encode_url = urlparse.urljoin(SERVER,f_obj.get('href'))
    logger.debug(encode_url)

    #stream=True avoids actually downloading the file, but it evaluates the redirection
    r = requests.get(encode_url, auth=(AUTHID,AUTHPW), headers={'content-type': 'application/json'}, allow_redirects=True, stream=True)
    try:
        r.raise_for_status
    except:
        logger.error('%s href does not resolve' %(f_obj.get('accession')))
        sys.exit()
    logger.debug(r)

    #this is the actual S3 https URL after redirection
    s3_url = r.url
    logger.debug(s3_url)

    #release the connection
    r.close()

    #split up the url into components
    o = urlparse.urlparse(s3_url)

    #pull out the filename
    filename = os.path.basename(o.path)

    #hack together the s3 cp url (with the s3 method instead of https)
    return filename, S3_SERVER.rstrip('/') + o.path


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
def noop(fastq):
    return {
        "report": None,
        "summary": None,
        "zip": None,
        "file": fastq
    }

@dxpy.entry_point("process")
def process(filename, bucket_url, project, folder):
    # Change the following to process whatever input this stage
    # receives.  You may also want to copy and paste the logic to download
    # and upload files here as well if this stage receives file input
    # and/or makes file output.

    logger.debug(filename)

    test = list( dxpy.find_data_objects(classname='file',
                           folder=folder, project=project, name_mode='exact',
                           name=filename, return_handler=False) )

    if not test or len(test) == 0:
        #cp the file from the bucket
        subprocess.check_call(shlex.split('aws s3 cp %s . --quiet' %(bucket_url)), stderr=subprocess.STDOUT)
        subprocess.check_call(shlex.split('ls -l %s' %(filename)))
        dx_file = dxpy.upload_local_file(filename, project=project, folder=folder)

    else:
        dx_file = test[0]
    reads_basename = filename.rstrip('.gz').rstrip('.fq').rstrip('.fastq')

    subprocess.check_call(['mkdir', 'output'])
    logger.info("Run QC")
    fqc_command = "/usr/bin/FastQC/fastqc " + filename + " -o output"
    logger.debug(fqc_command)
    stdio = subprocess.check_output(shlex.split(fqc_command))
    logger.debug(stdio)
    logger.debug(subprocess.check_output(['ls','-l', 'output']))
    subprocess.check_call(['unzip', "output/%s_fastqc.zip" % reads_basename])
    logger.info("Upload results")

    subprocess.check_call(['mv', 'fastq_fastqc/fastqc_data.txt', "%s_data.txt" % reads_basename])
    subprocess.check_call(['mv', 'fastq_fastqc/summary.txt', "%s_summary.txt" % reads_basename])

    report_dxfile = dxpy.upload_local_file("%s_data.txt" % reads_basename, folder=folder, project=project)
    summary_dxfile = dxpy.upload_local_file("%s_summary.txt" % reads_basename, folder=folder, project=project)
    zip_dxfile = dxpy.upload_local_file("output/%s_fastqc.zip" % reads_basename, folder=folder, project=project)
    logger.debug(report_dxfile)
    return {
        "file": dx_file,
        "report": report_dxfile,
        "summary": summary_dxfile,
        "zip": zip_dxfile
    }

@dxpy.entry_point("main")
def main(accession, key=None, debug=False):

    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.

    #files = [dxpy.DXFile(item) for item in files]

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    #for i, f in enumerate(files):
    #    dxpy.download_dxfile(f.get_id(), "files-" + str(i))

    # Split your work into parallel tasks.  As an example, the
    # following generates 10 subjobs running with the same dummy
    # input.

    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    (AUTHID,AUTHPW,SERVER) = processkey(key)

    url = SERVER + 'experiments/%s/?format=json&frame=embedded' %(accession)
    #get the experiment object
    logger.debug("%s - %s" % (url, AUTHID))
    response = encoded_get(url, AUTHID, AUTHPW)
    logger.debug(response)

    exp = response.json()
    # for some reason cannot write exp json to STDERR/logger

    '''
    Derive replicate structure and make directories
    '''
    project = dxpy.DXProject(os.environ['DX_PROJECT_CONTEXT_ID'])  ## should be default
    exp_folder = accession
    f = find_or_create_folder(project, exp_folder)
    for rep in exp['replicates']:
        rep_folder = "/%s_%s" % (rep['biological_replicate_number'], rep['technical_replicate_number'])
        rf = find_or_create_folder(project, rep_folder, root_folder='/'+exp_folder)


    subjobs = []
    for ff in exp['files']:
        if ff['file_format'] == 'fastq':
            folder = "/%s/%s_%s" % (exp_folder,
                ff['replicate']['biological_replicate_number'],
                ff['replicate']['technical_replicate_number'])
            file_name, bucket_url = get_bucket(SERVER, AUTHID, AUTHPW, ff)
            subjob_input = {
                "filename": file_name,
                "bucket_url": bucket_url,
                "project": project.get_id(),
                "folder": folder
            }
            subjobs.append(dxpy.new_dxjob(subjob_input, "process"))

    # The following line creates the job that will perform the
    # "postprocess" step of your app.  We've given it an input field
    # that is a list of job-based object references created from the
    # "process" jobs we just created.  Assuming those jobs have an
    # output field called "output", these values will be passed to the
    # "postprocess" job.  Because these values are not ready until the
    # "process" jobs finish, the "postprocess" job WILL NOT RUN until
    # all job-based object references have been resolved (i.e. the
    # jobs they reference have finished running).
    #
    # If you do not plan to have the "process" jobs create output that
    # the "postprocess" job will require, then you can explicitly list
    # the dependencies to wait for those jobs to finish by setting the
    # "depends_on" field to the list of subjobs to wait for (it
    # accepts either dxpy handlers or string IDs in the list).  We've
    # included this parameter in the line below as well for
    # completeness, though it is unnecessary if you are providing
    # job-based object references in the input that refer to the same
    # set of jobs.
    # If you would like to include any of the output fields from the
    # postprocess_job as the output of your app, you should return it
    # here using a job-based object reference.  If the output field in
    # the postprocess function is called "answer", you can pass that
    # on here as follows:
    #
    #return { "FastQC_reports": [ dxpy.dxlink(item) for item in postprocess_job.get_output_ref("report") ]}
    #
    # Tip: you can include in your output at this point any open
    # objects (such as gtables) which will be closed by a job that
    # finishes later.  The system will check to make sure that the
    # output object is closed and will attempt to clone it out as
    # output into the parent container only after all subjobs have
    # finished.

    output = {
                "files": [subjob.get_output_ref("file") for subjob in subjobs],
                "reports": [subjob.get_output_ref("report") for subjob in subjobs],
                "summaries": [subjob.get_output_ref("summary") for subjob in subjobs],
                "zips": [subjob.get_output_ref("zip") for subjob in subjobs],
    }


    return output

dxpy.run()
