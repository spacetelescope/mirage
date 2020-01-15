#! /usr/bin/env python
# coding=UTF-8
# Client for launching simulations for the Besancon Galaxy Model web service and Gravpot web service.
# Version 1.0,
# Date: 2019 january 14th.
# Author: Raphael Melior
# Institut UTINAM, OSU THETA, France
#
# .Change log:
#   - 2019-01-19:
#        model URL turned from http to https
#   - 2019-01-21
#        compatible python3 (add parentesis of print functions, long() => int())

# http://docs.python-requests.org/en/latest/user/quickstart/
# https://docs.python.org/2/library/xml.etree.elementtree.html

#________________________
#Example of usage for gravpot web service:
#   How to launch a simulation :
#./galmod_client.py --url "https://gravpot.utinam.cnrs.fr/" --user username --create -p RA 155.6 -p pmRA 2.5 -p Dec -12.5 -p pmDec -5.6 -p d 1.1 -p RV 15 --run
#
#   To see the status of a job :
#./galmod_client.py --url "https://gravpot.utinam.cnrs.fr/" --user username --job jobnumber
#
# see below (in paragraphe "main") various possible options and actions)
#________________________

import argparse
import datetime
import dateutil.parser
import requests # installs with : pip install requests
import xml.etree.ElementTree as ElementTree
import os
import getpass


XML_NS = {'uws':'http://www.ivoa.net/xml/UWS/v1.0', 'xlink':'http://www.w3.org/1999/xlink'}

class UWS:

    def __init__(self, url, auth):
        self._url = url
        self._auth = auth

    def getJobsList(self):
        r = requests.get(self._url+"/jobs", auth=self._auth)
        return self._getJobInfosFromReqResult(r, "get list of jobs")

    def createJob(self):
        r = requests.post(self._url+"/jobs", auth=self._auth)
        return self._getJobInfosFromReqResult(r, "creating job")

    def getJobInfos(self, idJob):
        r = requests.get(self._url+"/jobs/"+str(idJob), auth=self._auth)
        return self._getJobInfosFromReqResult(r, "get details of a job")

    def deleteJob(self, idJob):
        r = requests.delete(self._url+"/jobs/"+str(idJob), auth=self._auth)
        return self._getJobInfosFromReqResult(r, "deleting a job")

    def runJob(self, idJob):
        r = requests.post(self._url+"/jobs/"+str(idJob)+'/phase', auth=self._auth, data={'PHASE':'RUN'})
        return self._getJobInfosFromReqResult(r, "starting a job")

    def abortJob(self, idJob):
        r = requests.post(self._url+"/jobs/"+str(idJob)+'/phase', auth=self._auth, data={'PHASE':'ABORT'})
        return self._getJobInfosFromReqResult(r, "aborting a job")

    def maxDurationJob(self, idJob, seconds):
        r = requests.post(self._url+"/jobs/"+str(idJob)+'/executionduration', auth=self._auth, data={'EXECUTIONDURATION':seconds})
        return self._getJobInfosFromReqResult(r, "setting execution duration of a job")

    def timeDestructJob(self, idJob, timeStr):
        r = requests.post(self._url+"/jobs/"+str(idJob)+'/destruction', auth=self._auth, data={'DESTRUCTION':timeStr})
        return self._getJobInfosFromReqResult(r, "setting execution duration of a job")

    def setParams(self, idJob, params):
        #if len(params) == 1:   # TODO
        #   r = requests.put(self._url+"/jobs/"+str(idJob)+'/parameters/'+params[0][0], data=params[0][1])
        #   print(r.text)
        #   return getJobInfos(idJob)#self._getJobInfosFromReqResult(r, "setting parameters of a job")
        #elif len(params) > 1:
            data={}
            for p in params:
                data[p[0]] = p[1]
            r = requests.post(self._url+"/jobs/"+str(idJob)+'/parameters', auth=self._auth, data=data)
            return self._getJobInfosFromReqResult(r, "setting parameters of a job")

    def _getJobInfosFromReqResult(self, r, phase):
        if (len(r.history) > 0):
            print("SERVER RESPONSE : "+r.history[0].text)

        if (r.status_code == 200):
            #et = ElementTree
            #ElementTree.register_namespace('uws', 'http://www.ivoa.net/xml/UWS/v1.0')
            for ns in XML_NS:
                ElementTree.register_namespace(ns, XML_NS[ns])

            try:
                return ElementTree.fromstring(r.text)
            except Exception as e:
                print("Invalid response :")
                print(r.text)
                raise e
        else:
            raise JobError(r, phase)




# Convertion into text

def xmlTextIfFound(xmlElem):
    if (xmlElem == None) or (xmlElem.text == None):
        return ''
    else:
        return xmlElem.text

def xmlDateIfFound(xmlElem):
    if (xmlElem == None) or (xmlElem.text == None):
        return '-'
    else:
        d = dateutil.parser.parse(xmlElem.text).astimezone(dateutil.tz.tzutc())
        return d.strftime("%Y-%m-%d %H:%M")


# Exceptions

class JobError(Exception):
    def __init__(self, http, when):
        self.status = http.status_code
        self.explain = http.text
        self.when = when
    def __str__(self):
        return "Error when "+self.when+" : "+str(self.status)+" "+self.explain



# Main

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Makes requests to web services, either the Besancon Model of the Galaxy or Gravpot web service', add_help=False)

    serverGrp = parser.add_argument_group('Server options', 'Options for server connection.')
    serverGrp.add_argument('--url', help='URL of the server.', default='https://model.obs-besancon.fr/ws')
    serverGrp.add_argument('--user', help='Username used to connect to server, by default it is the system username ('+os.getlogin()+').', default=os.getlogin())
    serverGrp.add_argument('--pass', dest='passwd', help='Password used to connect to server, by default the password is prompted. Usage of this option is discouraged because it will be visible, for example, with "ps" command and written to history.')

    jobGrp = parser.add_argument_group('Job selection', 'Select the job you want to manipulate. To create a new one type "--create". If this option is not present, it returns your list of jobs and exit ignoring options in "Actions".').add_mutually_exclusive_group()
    jobGrp.add_argument('-c', '--create', help='Create a new job.', action='store_true')
    jobGrp.add_argument('-j', '--job', type=int, help='The job affected by this request.')

    actionGrp = parser.add_argument_group('Actions', 'Configure and make actions on the selected job.')
    actionGrp.add_argument('-p', '--param', help='Set parameter named PARAMNAME with VALUE.', action='append', nargs=2, metavar=('PARAMNAME', 'VALUE'))
    actionGrp.add_argument('--execdur', help='Set maximum execution duration of a job in seconds (this setting can be overridden by server\'s configuration).', type=int)
    actionGrp.add_argument('--tdest', help='Set destruction time of a job (date+time in ISO8601 format, this setting can be overridden by server\'s configuration).')
    actionGrp.add_argument('--run', help='Send job for computation.', action='store_true')
    actionGrp.add_argument('--abort', help='Abort the job.', action='store_true')
    actionGrp.add_argument('--delete', help='Delete the job.', action='store_true')

    miscGrp = parser.add_argument_group('Miscellaneous')
    miscGrp.add_argument('-h', '--help', help='Show this help message and exit.', action='store_true')

    args = parser.parse_args()

    if (args.help):
        parser.print_help()
        exit(0)

    if (args.passwd != None):
        passwd = args.passwd
    else:
        passwd = getpass.getpass('Password for '+args.user+': ')

    #print(args.user+' '+passwd)
    uws = UWS(args.url, (args.user, passwd))

    #print(args) #args.output
    job = None
    if (args.create):
        try:
            print("Creating job")
            job = uws.createJob()
        except JobError as e:
            print(e)
            exit(2)

    elif (args.job != None):
        try:
            job = uws.getJobInfos(args.job)
        except JobError as e:
            print(e)
            exit(2)

    else:
        print(" JobID        phase       startTime          endTime            destruction    ")
        print("----------   ---------   ----------------   ----------------   ----------------")
        #print("0123456789   COMPLETED   2009-05-19 17:12   2009-05-19 17:15   2009-06-06 17:12")
        jobs = uws.getJobsList()
        for job in jobs.findall('uws:job', XML_NS):
            jobId = xmlTextIfFound(job.find('uws:jobId', XML_NS))
            phase = xmlTextIfFound(job.find('uws:phase', XML_NS))
            startTime = xmlDateIfFound(job.find('uws:startTime', XML_NS))
            endTime = xmlDateIfFound(job.find('uws:endTime', XML_NS))
            destruction = xmlDateIfFound(job.find('uws:destruction', XML_NS))
            print("%10d   %9s   %16s   %16s   %16s" % (int(jobId), phase, startTime, endTime, destruction))
        exit(0)


    idJob = int(job.find('uws:jobId', XML_NS).text)

    if args.param:
        job = uws.setParams(idJob, args.param)

    if args.execdur:
        job = uws.maxDurationJob(idJob, args.execdur)

    if args.tdest:
        job = uws.timeDestructJob(idJob, args.tdest)

    if args.run:
        job = uws.runJob(idJob)
    elif args.abort:
        job = uws.abortJob(idJob)

    print("id                :"+str(idJob))
    print("owner             :"+xmlTextIfFound(job.find('uws:ownerId', XML_NS)))
    print("phase             :"+xmlTextIfFound(job.find('uws:phase', XML_NS)))
    print("startTime         :"+xmlTextIfFound(job.find('uws:startTime', XML_NS)))
    print("endTime           :"+xmlTextIfFound(job.find('uws:endTime', XML_NS)))
    print("destruction       :"+xmlTextIfFound(job.find('uws:destruction', XML_NS)))
    print("executionDuration :"+xmlTextIfFound(job.find('uws:executionDuration', XML_NS)))

    parameters = job.find('uws:parameters', XML_NS)
    if parameters is not None:
        print('\nParameters ........')
        for param in parameters.findall('uws:parameter', XML_NS):
            print('{:<18}:{}'.format(param.get('id'), xmlTextIfFound(param)))

    errorSummary = job.find('uws:errorSummary', XML_NS)
    if errorSummary is not None:
        print('\nErrors ............')
        for err in errorSummary.findall('uws:message', XML_NS):
            print('error             :'+xmlTextIfFound(err))

    results = job.find('uws:results', XML_NS)
    if results is not None:
        print('\nResults ...........')
        for res in results.findall('uws:result', XML_NS):
            print('{:<18}:{}'.format(res.get('id'), res.get('{http://www.w3.org/1999/xlink}href')))

    if args.delete:
        uws.deleteJob(idJob)
