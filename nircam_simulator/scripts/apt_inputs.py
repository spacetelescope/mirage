#! /usr/bin/env python

'''
Given APT output files, read in data relevant to the data simulator,
organize, and create input files for the simulator.

Inputs:

xml file - Name of xml file exported from APT.
pointing file - Name of associated pointing file exported from APT.
siaf - Name of csv version of SIAF.

Optional inputs:

epoch_list - Name of ascii file which lists observation labels and 
             associated starting observation time, as well as telescope
             roll angle (PAV3). If you wish to have observations back-
             to-back, give them all the same starting time.

Outputs:

output_csv - Ascii table containing all relevant information needed
             to construct a ramp simulator input file for each
             observation/dither/detector combination.

Dependencies:

argparse, lxml, astropy, numpy, collections

JWST Calibration pipeline (only for the use of set_telescope_pointing.py)
in order to translate PAV3 values to local roll angles for each detector.

rotations.py - Colin Cox's module of WCS-related functions


HISTORY:

July 2017 - V0: Initial version. Bryan Hilbert

'''
import os
import sys
import collections
import argparse
from lxml import etree
from astropy.table import Table, Column
from astropy.io import ascii
import numpy as np
import yaml
from . import rotations
from . import set_telescope_pointing_separated as set_telescope_pointing

class AptInput:
    def __init__(self):
        self.input_xml = '' #e.g. 'GOODSS_ditheredDatasetTest.xml'
        self.output_csv = None #e.g. 'GOODSS_ditheredDatasetTest.csv'
        self.pointing_file = '' #e.g. 'GOODSS_ditheredDatasetTest.pointing'
        self.siaf = ''
        self.observation_table = ''
        self.output_csv = None
        

    def read_xml(self,infile):
        #read in APT xml file. WFSS and/or direct exposures can be
        #in the file

        #open file, get tree
        with open(infile) as f:
            tree = etree.parse(f)
        namespaces = tree.getroot().nsmap.copy()
        namespaces['apt'] = namespaces[None]
        del namespaces[None]

        #creat dictionary to hold all relevant exposure info
        dict = {}
        dict['PI_Name'] = []
        dict['ProposalID'] = []
        dict['Title'] = []
        dict['Proposal_category'] = []
        dict['Science_category'] = []
        dict['Mode'] = []
        dict['Module'] = []
        dict['Subarray'] = []
        dict['PrimaryDitherType'] = []
        dict['PrimaryDithers'] = []
        dict['SubpixelDitherType'] = []
        dict['SubpixelPositions'] = []
        dict['ShortFilter'] = []
        dict['LongFilter'] = []
        dict['ReadoutPattern'] = []
        dict['Groups'] = []
        dict['Integrations'] = []
        dict['ShortPupil'] = []
        dict['LongPupil'] = []
        dict['Grism'] = []
        dict['CoordinatedParallel'] = []

        #get high level information: proposal info
        #P.I. Name, etc
        propid_default = 42424
        proptitle_default = 'Looking for my towel'
        scicat_default = 'Planets and Planet Formation'
        piname_default = 'D.N. Adams'

        apt = '{http://www.stsci.edu/JWST/APT}'
        propinfo = tree.find(apt+'ProposalInformation')
        try:
            proptitle = propinfo.find(apt+'Title').text
        except:
            proptitle = proptitle_default
        try:
            propid = propinfo.find(apt+'ProposalID').text
        except:
            propid = propid_default
        propcat = 'GO'
        try:
            scicat = propinfo.find(apt+'ScientificCategory').text
        except:
            scicat = scicat_default
        try:
            piinfo = propinfo.find(apt+'PrincipalInvestigator')
            piadd = piinfo.find(apt+'InvestigatorAddress')
            firstname = piadd.find(apt+'FirstName').text
            lastname = piadd.find(apt+'LastName').text
            piname = firstname + ' ' + lastname
        except:
            piname = piname_default


        #Also look at mosaic details. All we really need to know
        #is how many tiles will be taken.
        
        
        #more streamlined version with valid (i belive) assumptions
        #about format of xml file:
        #within an observation: only 1 instrument, 1 template, 1 nci:NircamImaging or 1 ncwfss:NircamWfss
        obspath = "//apt:Observation"
        obsresults = tree.xpath(obspath, namespaces=namespaces)
        obsgrouppath = "//apt:ObservationGroup"
        obsgroupresults = tree.xpath(obsgrouppath,namespaces=namespaces)

        apt = '{http://www.stsci.edu/JWST/APT}'
        nci = "{http://www.stsci.edu/JWST/APT/Template/NircamImaging}"
        ncwfss = "{http://www.stsci.edu/JWST/APT/Template/NircamWfss}"
        mos = "{http://www.stsci.edu/JWST/APT/MosaicParameters}"
        for obs in obsresults:
            obs_tuple_list = []
            numbele = obs.find(apt+'Number')
            number = numbele.text
            labelele = obs.find(apt+'Label')
            if labelele is not None:
                label = labelele.text
            else:
                label = 'None'
            targele = obs.find(apt+'TargetID')
            targsplit = targele.text.split()
            target = ''
            for i in range(1,len(targsplit)):
                target += targsplit[i]
            instele = obs.find(apt+'Instrument')

            if instele.text == 'NIRCAM':
                template = obs.find(apt+'Template')
                coordparallel = obs.find(apt+'CoordinatedParallel').text
                imaging_temp = template.find(nci+'NircamImaging')
                if imaging_temp is not None:
                    typeflag = 'Imaging'
                    grismval = 'N/A'
                    short_pupil = 'CLEAR'
                    mod = imaging_temp.find(nci+'Module').text
                    subarr = imaging_temp.find(nci+'Subarray').text
                    pdithtype = imaging_temp.find(nci+'PrimaryDitherType').text
                    try:
                        pdither = imaging_temp.find(nci+'PrimaryDithers').text
                    except:
                        pdither = '1'
                    sdithtype = imaging_temp.find(nci+'SubpixelDitherType').text
                    try:
                        sdither = imaging_temp.find(nci+'SubpixelPositions').text
                    except:
                        try:
                            stemp = imaging_temp.find(nci+'CoordinatedParallelSubpixelPositions').text
                            sdither = np.int(stemp[0])
                        except:
                            sdither = '1'
                    filtele = imaging_temp.find(nci+'Filters') 
                    filtconfigeles = filtele.findall(nci+'FilterConfig')
                    for fcele in filtconfigeles:
                        sfilt = fcele.find(nci+'ShortFilter').text
                        lfilt = fcele.find(nci+'LongFilter').text
                        rpatt = fcele.find(nci+'ReadoutPattern').text
                        grps = fcele.find(nci+'Groups').text
                        ints = fcele.find(nci+'Integrations').text
                
                        #separate pupil and filter in case of filter that is 
                        #mounted in the pupil wheel
                        if '+' in sfilt:
                            p = sfilt.find('+')
                            short_pupil = sfilt[0:p]
                            sfilt = sfilt[p+1:]
                        else:
                            short_pupil = 'CLEAR'
                    
                        if '+' in lfilt:
                            p = lfilt.find('+')
                            long_pupil = lfilt[0:p]
                            lfilt = lfilt[p+1:]
                        else:
                            long_pupil = 'CLEAR'
                
                        tup_to_add = (piname,propid,proptitle,propcat,
                                      scicat,typeflag,mod,subarr,pdithtype,
                                      pdither,sdithtype,sdither,sfilt,lfilt,
                                      rpatt,grps,ints,short_pupil,
                                      long_pupil,grismval,coordparallel)
                        dict = self.add_exposure(dict,tup_to_add)
                        obs_tuple_list.append(tup_to_add)
                        
                wfss_temp = template.find(ncwfss+'NircamWfss')
                if wfss_temp is not None:
                    mod = wfss_temp.find(ncwfss+'Module').text
                    subarr = wfss_temp.find(ncwfss+'Subarray').text
                    grismval = wfss_temp.find(ncwfss+'Grism').text
                    if grismval == 'BOTH':
                        grismval = ['GRISMR','GRISMC']
                    else:
                        grismval = [grismval]
                    #pdithtype = wfss_temp.find(ncwfss+'PrimaryDitherType').text
                    #pdither = wfss_temp.find(ncwfss+'PrimaryDithers').text
                    #sdither = wfss_temp.find(ncwfss+'SubpixelPositions').text
                    #sdithtype = wfss_temp.find(ncwfss+'SubpixelPositions').text
                    explist = wfss_temp.find(ncwfss+'ExposureList')
                    expseqs = explist.findall(ncwfss+'ExposureSequences')
            
                    #if BOTH was specified for the grism,
                    #then we need to repeat the sequence of
                    #grism/direct/grism/direct/outoffield for each grism
                    for gnum in range(len(grismval)):
                        for expseq in expseqs:
                            #sequence = grism,direct,grism,direct,outoffield
                            #if grism == both, sequence is done for grismr,
                            #then repeated for grismc
                            grismvalue = grismval[gnum]
                            #need to switch the order of the grism and direct
                            #exposures in order for them to be chronological
                            grismexp = expseq.find(ncwfss+'GrismExposure')
                            typeflag = 'WFSS'
                            sfilt = grismexp.find(ncwfss+'ShortFilter').text
                            lfilt = grismexp.find(ncwfss+'LongFilter').text
                            rpatt = grismexp.find(ncwfss+'ReadoutPattern').text
                            groups = grismexp.find(ncwfss+'Groups').text
                            integrations = grismexp.find(ncwfss+'Integrations').text

                            pdithtype = wfss_temp.find(ncwfss+'PrimaryDitherType').text
                            pdither = wfss_temp.find(ncwfss+'PrimaryDithers').text
                            sdither = wfss_temp.find(ncwfss+'SubpixelPositions').text
                            sdithtype = wfss_temp.find(ncwfss+'SubpixelPositions').text

                            
                            #separate pupil and filter in case of filter
                            #that is mounted in the pupil wheel
                            if '+' in sfilt:
                                p = sfilt.find('+')
                                short_pupil = sfilt[0:p]
                                sfilt = sfilt[p+1:]
                            else:
                                short_pupil = 'CLEAR'
                    
                            long_pupil = grismvalue
                            tup_to_add = (piname,propid,proptitle,propcat,
                                          scicat,typeflag,mod,subarr,
                                          pdithtype,pdither,sdithtype,
                                          sdither,sfilt,lfilt,rpatt,groups,
                                          integrations,short_pupil,long_pupil,
                                          grismvalue,coordparallel)

                            dict = self.add_exposure(dict,tup_to_add)
                            obs_tuple_list.append(tup_to_add)
                            
                            directexp = expseq.find(ncwfss+'DiExposure')
                            typeflag = 'Imaging'
                            pdither = '1' #direct image has no dithers
                            sdither = '1' #direct image has no dithers
                            sdithtype = '1' #direct image has no dithers
                            grismvalue = 'N/A'
                            sfilt = directexp.find(ncwfss+'ShortFilter').text
                            lfilt = directexp.find(ncwfss+'LongFilter').text
                            rpatt = directexp.find(ncwfss+'ReadoutPattern').text
                            grps = directexp.find(ncwfss+'Groups').text
                            ints = directexp.find(ncwfss+'Integrations').text
                
                            #separate pupil and filter in case of filter
                            #that is mounted in the pupil wheel
                            if '+' in sfilt:
                                p = sfilt.find('+')
                                short_pupil = sfilt[0:p]
                                sfilt = sfilt[p+1:]
                            else:
                                short_pupil = 'CLEAR'
                    
                            if '+' in lfilt:
                                p = lfilt.find('+')
                                long_pupil = lfilt[0:p]
                                lfilt = lfilt[p+1:]
                            else:
                                long_pupil = 'CLEAR'
            
                            direct_tup_to_add = (piname,propid,proptitle,propcat,
                                                 scicat,typeflag,mod,subarr,pdithtype,
                                                 pdither,sdithtype,sdither,sfilt,lfilt,
                                                 rpatt,grps,ints,short_pupil,long_pupil,
                                                 grismvalue,coordparallel)
                            dict = self.add_exposure(dict,direct_tup_to_add)
                            obs_tuple_list.append(tup_to_add)
                            
                        #Now we need to add the two out-of-field exposures, which are
                        #not present in the APT file (but are in the associated pointing
                        #file from APT. We can just
                        #duplicate the entries for the direct images taken immediately
                        #prior. BUT, will there ever be a case where there is no preceding
                        #direct image?
                        dict = self.add_exposure(dict,direct_tup_to_add)
                        dict = self.add_exposure(dict,direct_tup_to_add)
                        obs_tuple_list.append(tup_to_add)
                        obs_tuple_list.append(tup_to_add)

                #Now we need to look for mosaic details, if any
                mosaicele = obs.find(apt+'MosaicParameters')
                mostile = mosaicele.findall(apt+'MosaicTiles')
                numtiles = len(mostile)
                if numtiles > 1:
                    print("Found {} mosaic tiles for observation {}".format(numtiles,obs))
                    for i in range(numtiles-1):
                        for tup in obs_tuple_list:
                            dict = self.add_exposure(dict,tup)
        return dict

            
    def add_exposure(self,dictionary,tup):
        #add an exposure to the dictionary 
        dictionary['PI_Name'].append(tup[0])
        dictionary['ProposalID'].append(tup[1])
        dictionary['Title'].append(tup[2])
        dictionary['Proposal_category'].append(tup[3])
        dictionary['Science_category'].append(tup[4])
        dictionary['Mode'].append(tup[5])
        dictionary['Module'].append(tup[6])
        dictionary['Subarray'].append(tup[7])
        dictionary['PrimaryDitherType'].append(tup[8])
        dictionary['PrimaryDithers'].append(tup[9])
        dictionary['SubpixelDitherType'].append(tup[10])
        dictionary['SubpixelPositions'].append(tup[11])
        dictionary['ShortFilter'].append(tup[12])
        dictionary['LongFilter'].append(tup[13])
        dictionary['ReadoutPattern'].append(tup[14])
        dictionary['Groups'].append(tup[15])
        dictionary['Integrations'].append(tup[16])
        dictionary['ShortPupil'].append(tup[17])
        dictionary['LongPupil'].append(tup[18])
        dictionary['Grism'].append(tup[19])
        dictionary['CoordinatedParallel'].append(tup[20])
        return dictionary
        
    def read_wfss_xml(self,infile):
        #read APT xml file for WFSS mode observations
        #first, set up variables 
        MyList = collections.OrderedDict()
        MyList['Module'] = []
        MyList['Subarray'] = [] 
        MyList['Grism'] = []
        MyList['PrimaryDitherType'] = [] 
        MyList['PrimaryDithers'] = []
        MyList['SubpixelPositions'] = []
        MyList['TargID'] = []
        MyFilterList = collections.OrderedDict()
        MyFilterList['Module'] = []
        MyFilterList['Subarray'] = [] 
        MyFilterList['Grism'] = []
        MyFilterList['PrimaryDitherType'] = [] 
        MyFilterList['PrimaryDithers'] = []
        MyFilterList['SubpixelPositions'] = []
        MyFilterList['Mode'] = []
        MyFilterList['ShortFilter'] = []
        MyFilterList['LongFilter'] = []
        MyFilterList['ReadoutPattern'] = []
        MyFilterList['Groups'] = []
        MyFilterList['Integrations'] = []
        MyTargList = []

        #read in the full file
        f = open(self.input_xml)
        fullfile = f.readlines()
        f.close()
        
        #now find the lines corresponding to the beginning of each
        #exposure list.
        wfss_start = np.array([]).astype(np.int)
        wfss_end = np.array([]).astype(np.int)
        explist_start = np.array([]).astype(np.int)
        directlist_start = np.array([]).astype(np.int)
        grismlist_start = np.array([]).astype(np.int)
        targlines = np.array([]).astype(np.int)

        #default values in case of missing data in APT file
        propid = '42424' 
        title = 'I need to find my towel'
        piname = 'D.N. Adams'
        pistart = 0
        piend = -1
        prop_category = 'GO'
        science_category = 'extrasolar towels'
        for linenum in range(len(fullfile)):
            if "<ncwfss:NircamWfss>" in fullfile[linenum]:
                wfss_start = np.append(wfss_start,linenum)
            if "</ncwfss:NircamWfss>" in fullfile[linenum]:
                wfss_end = np.append(wfss_end,linenum)
            if "<ncwfss:ExposureList>" in fullfile[linenum]:
                explist_start = np.append(explist_start,linenum)
            if "<ncwfss:DiExposure>" in fullfile[linenum]:
                directlist_start = np.append(directlist_start,linenum)
            if "<ncwfss:GrismExposure>" in fullfile[linenum]:
                grismlist_start = np.append(grismlist_start,linenum)
            #get the proposal ID number
            if "<ProposalID>" in fullfile[linenum]:
                propid = self.extract_value(fullfile[linenum])
            if "<Title>" in fullfile[linenum]:
                title = self.extract_value(fullfile[linenum])
            if "<PrincipalInvestigator>" in fullfile[linenum]:
                pistart = linenum
            if "</PrincipalInvestigator>" in fullfile[linenum]:
                piend = linenum
            if "<ProposalCategory>" in fullfile[linenum]:
                prop_category = self.extract_value(fullfile[linenum+1])[:-1]
            if "<ScientificCategory>" in fullfile[linenum]:
                science_category = self.extract_value(fullfile[linenum])  
            #if "<Target ID>" in fullfile[linenum]:
            #    targlines.append(linenum)

        if pistart > 0:
            for lnum in range(pistart,piend):
                if "<FirstName>" in fullfile[lnum]:
                    first = self.extract_value(fullfile[lnum])
                if "<LastName>" in fullfile[lnum]:
                    last = self.extract_value(fullfile[lnum])
            piname = first + ' ' + last
                
        #now, work on each wfss_start entry individually.
        #each one of these will have grism exposures, optional
        #direct exposure (singular), and out of field exposures (2).
        #We need to keep these all grouped together so that we end
        #with an exposure list that is in chronological order
        for wfssstart,wfssend in zip(wfss_start,wfss_end):
            for listele,addline in zip(MyList.keys(),range(1,7)):
                gt = fullfile[wfssstart+addline].find('>')
                lt = fullfile[wfssstart+addline].find('<',gt)
                MyList[listele] = fullfile[wfssstart+addline][gt+1:lt]
            #associate a target ID with each
            #prevtarg = np.where(targlines < wfssstart)
            #targ = self.extract_value(fullfile[prevtarg][-1]).split()
            #fulltarg = ''
            #for i in range(1,len(targ)):
            #    fulltarg = fulltarg + targ[i]
            ##MyTargList.append(fulltarg.strip())

            #now get info on the grism and optional direct images
            #that are only within the current wfss_start entry
    
            #first get all the grism exposure info
            gline = ((grismlist_start > wfssstart) &
                     (grismlist_start < wfssend))

            #make sure there is a WFSS entry in this exposure list
            if np.sum(gline) > 0:
                gentries = grismlist_start[gline]
        
                #loop over grism entry start lines
                for gindex in range(len(gentries)):
                    grismstart = gentries[gindex]
                    if gindex != 0:
                        prev_grism = gentries[gindex-1]
                    else:
                        prev_grism = wfssstart
                
                    MyFilterList['Mode'].append('WFSS')
                    for listele,addline in zip(MyFilterList.keys()[7:12],range(7,12)):
                        gt = fullfile[grismstart+addline-6].find('>')
                        lt = fullfile[grismstart+addline-6].find('<',gt)
                        MyFilterList[listele].append(fullfile[grismstart+addline-6][gt+1:lt])
                    for key in MyList:
                        MyFilterList[key].append(MyList[key])
                    MyTargList.append(fulltarg.strip())
                    #now get any direct exposure info that is tied to this 
                    #grism exposure
                    dline = ((directlist_start > prev_grism) &
                             (directlist_start < grismstart))
            
                    #make sure there is a direct image entry in this exposure list
                    if np.sum(dline) > 0:
                        for directstart in directlist_start[dline]:
                            MyFilterList['Mode'].append('Imaging')
                            MyTargList.append(fulltarg.strip())
                            for listele,addline in zip(MyFilterList.keys()[7:12],range(7,12)):
                                gt = fullfile[directstart+addline-6].find('>')
                                lt = fullfile[directstart+addline-6].find('<',gt)
                                MyFilterList[listele].append(fullfile[directstart+addline-6][gt+1:lt])
                            for key in MyList:
                                if key not in ['PrimaryDithers','SubpixelPositions']:
                                    MyFilterList[key].append(MyList[key])
                                else:
                                    MyFilterList[key].append('1')
            #now we need to add the two OUT OF FIELD expoures,
            #which are not in the xml file. They use the same
            #readout pattern/groups/ints as the direct image.
            #This is seen within APT itself, but is not in the
            #xml file. Since everything is the same as the direct
            #image taken immediately prior, we can just duplicate
            #the dictionary entries for the direct image.
            for key in MyFilterList:
                MyFilterList[key].extend((MyFilterList[key][-1],MyFilterList[key][-1]))

        #add proposal info
        MyFilterList['ProposalID'] = []
        if len(MyFilterList['Module']) > 0:
            MyFilterList['ProposalID'] = [np.int(propid)]*len(MyFilterList['Module'])
            MyFilterList['Title'] = [title]*len(MyFilterList['Module'])
            MyFilterList['PI_Name'] = [piname]*len(MyFilterList['Module'])
            MyFilterList['Proposal_category'] = [prop_category]*len(MyFilterList['Module'])
            MyFilterList['Science_category'] = [science_category]*len(MyFilterList['Module'])
            
        #now we need to deal with the pupil values.
        swpupillist = ['CLEAR'] * len(MyFilterList['Mode'])
        lwpupillist = ['CLEAR'] * len(MyFilterList['Mode'])
        for i in range(len(MyFilterList['Mode'])):
            if MyFilterList['Mode'][i] == 'WFSS':
                lwpupillist[i] = MyFilterList['Grism'][i]
            else:
                if '+' in MyFilterList['LongFilter'][i]:
                    p = MyFilterList['LongFilter'][i].find('+')
                    pup = MyFitlerList['LongFilter'][i][0:p]
                    f1 = MyFitlerList['LongFilter'][i][p+1:]
                    lwpupillist[i] = pup
                    MyFilterList['LongFilter'][i] = f1
            if '+' in MyFilterList['ShortFilter'][i]:
                p = MyFilterList['ShortFilter'][i].find('+')
                pup = MyFilterList['ShortFilter'][i][0:p]
                f1 = MyFilterList['ShortFilter'][i][p+1:]
                swpupillist[i] = pup
                MyFilterList['ShortFilter'][i] = f1
        MyFilterList['ShortPupil'] = swpupillist
        MyFilterList['LongPupil'] = lwpupillist

        #add in target names
        #print('list lengths: ',len(MyTargList),len(MyFilterList['ShortPupil'])))
        #stop
        #MyFilterList['TargID'] = MyTargList

        #add subpixeldithertype, to be consistent with imaging output
        MyFilterList['SubpixelDitherType'] = MyFilterList['SubpixelPositions']
        return MyFilterList


    def extract_value(self,line):
        #extract text from xml line
        gt = line.find('>')
        lt = line.find('<',gt)
        return line[gt+1:lt]
                                                  
        
    def read_imaging_xml(self,infile):
        #read APT xml file for imaging mode obs

        #first, a cheat. get proposal id by reading in file as ascii,
        #because I can't figure out the xml way to do it
                #read in the full file
        f = open(self.input_xml)
        fullfile = f.readlines()
        f.close()
        
        #get proposal information
        #default values in case of missing data in APT file
        propid = '42424' 
        title = 'Looking for my towel'
        piname = 'D.N. Adams'
        pistart = 0
        piend = -1
        prop_category = 'GO'
        science_category = 'extrasolar towels'

        for linenum in range(len(fullfile)):
            if "<ProposalID>" in fullfile[linenum]:
                propid = self.extract_value(fullfile[linenum])
            if "<Title>" in fullfile[linenum]:
                title = self.extract_value(fullfile[linenum])
            if "<PrincipalInvestigator>" in fullfile[linenum]:
                pistart = linenum
            if "</PrincipalInvestigator>" in fullfile[linenum]:
                piend = linenum
            if "<ProposalCategory>" in fullfile[linenum]:
                lt = fullfile[linenum+1].find('<')
                gt = fullfile[linenum+1].find('>')
                prop_category = fullfile[linenum+1][lt+1:gt][:-1]
            if "<ScientificCategory>" in fullfile[linenum]:
                science_category = self.extract_value(fullfile[linenum])             

        if pistart > 0:
            for lnum in range(pistart,piend):
                if "<FirstName>" in fullfile[lnum]:
                    first = self.extract_value(fullfile[lnum])
                if "<LastName>" in fullfile[lnum]:
                    last = self.extract_value(fullfile[lnum])
            piname = first + ' ' + last

                
        path = "//apt:Observation[apt:Instrument[contains(string(), '{}')]]/apt:Template/nci:NircamImaging".format('NIRCAM')
        
        targpath = "//apt:Observation"

        # READ XML file
        with open(infile) as f:
            tree = etree.parse(f)

        # APT makes extensive use of XML namespaces
        # (e.g. 'xmlns:nsmsasd="http://www.stsci.edu/JWST/APT/Template/NirspecMSAShortDetect"')
        # so we have to as well
        namespaces = tree.getroot().nsmap.copy()
        # There is no 'default' namespace for XPath (used below), but the lxml parser
        # does respect a default namespace, so we have to update its name from
        # 'None' to 'apt'
        namespaces['apt'] = namespaces[None]
        del namespaces[None]

        # Find your specific Observation
        results = tree.xpath(path, namespaces=namespaces)
        targresults = tree.xpath(targpath, namespaces=namespaces)

        #set up variables for output
        MyList = {'Module': [], 'Subarray': [], 'PrimaryDitherType': [],
                  'PrimaryDithers': [], 'SubpixelDitherType': [],
                  'SubpixelPositions': []}
        MyFilterList = {'ShortFilter': [], 'LongFilter': [],
                        'ReadoutPattern': [], 'Groups': [], 'Integrations': []}
        finalList = {'Mode':[], 'Module': [], 'Subarray': [],
                     'PrimaryDitherType': [], 'PrimaryDithers': [],
                     'SubpixelDitherType': [], 'SubpixelPositions': [],
                     'ShortFilter': [], 'LongFilter': [], 'ReadoutPattern': [],
                     'Groups': [], 'Integrations': []}
        
        for ExposureList in results:
            #reset the lists in the dictionaries to be empty at
            #the beginning of each Template
            MyList = {'Module': [], 'Subarray': [], 'PrimaryDitherType': [],
                      'PrimaryDithers': [], 'SubpixelDitherType': [],
                      'SubpixelPositions': []}
            MyFilterList = {'ShortFilter': [], 'LongFilter': [],
                            'ReadoutPattern': [], 'Groups': [], 'Integrations': []}

            for item in MyList:
                entryList = ExposureList.xpath('nci:%s' % item,namespaces=namespaces)
                for entry in entryList:
                    MyList[item].append(entry.text)
            for item in MyFilterList:
                entryList = ExposureList.xpath('nci:Filters/nci:FilterConfig/nci:%s' % item,namespaces=namespaces)
                for entry in entryList:
                    MyFilterList[item].append(entry.text)
            
            #Add in a mode keyword so that we can easily separate
            #imaging from wfss entries. This will be useful once
            #these outputs are passed to the tool for making simulator
            #input files
            MyList['Mode'] = ['Imaging']*len(MyList['Module'])
                       
            #duplicate entries in the MyList dictionary so that the length
            #matches the myFilterList dictionary
            n_module = len(MyList['Module'])
            n_filter = len(MyFilterList['ShortFilter'])
            reps = n_filter - n_module
            for key in MyFilterList:
                finalList[key] = finalList[key] + MyFilterList[key]
            for key in MyList:
                finalList[key] = finalList[key] + MyList[key]*(reps+1)

        #check the filters. In the case where a pupil wheel-mounted filter
        #is used, the filter name will be "filter1+filter2". Separate into
        #filter and pupil entries
        shortplist = ['CLEAR']*len(finalList['ShortFilter'])
        longplist = ['CLEAR']*len(finalList['LongFilter'])
        for key in ['ShortFilter','LongFilter']:
            for i in range(len(finalList['ShortFilter'])):
                filt = finalList[key][i]
                if '+' in filt:
                    p = filt.find('+')
                    pupil = filt[0:p]
                    f1 = filt[p+1:]
                    if key == 'ShortFilter':
                        finalList[key][i] = f1                
                        shortplist[i] = pupil
                    if key == 'LongFilter':
                        finalList[key][i] = f1
                        longplist[i] = pupil
        finalList['ShortPupil'] = shortplist
        finalList['LongPupil'] = longplist

        #for consistency with the output from the WFSS reader
        finalList['Grism'] = ['N/A'] * len(finalList['Mode'])

        #add proposal info lines
        finalList['ProposalID'] = [np.int(propid)]*len(finalList['Module'])
        finalList['Title'] = [title]*len(finalList['Module'])
        finalList['PI_Name'] =  [piname]*len(finalList['Module'])
        finalList['Proposal_category'] = [prop_category]*len(finalList['Module'])
        finalList['Science_category'] = [science_category]*len(finalList['Module'])


        #for code checking
        ascii.write(Table(finalList), 'test_imaging.csv', format='csv', overwrite=True)
        
        return finalList


    def expand_for_dithers(self,indict):
        #Expand a given dictionary to create one entry
        #for each dither
        #define the dictionary to hold the expanded entries

        #in here we should also reset the primary and subpixel dither
        #numbers to 1, to avoid confusion.
        
        expanded = {}
        for key in indict:
            expanded[key] = []
    
        #loop over entries in dict and duplicate by the
        #number of dither positions  
        #keys = np.array(indict.keys())
        keys = indict.keys()
        for i in range(len(indict['PrimaryDithers'])):
            #entry = np.array([item[i] for item in dict.values()])
            arr = np.array([item[i] for item in indict.values()])
            entry = dict(zip(keys,arr))
            
            #subpix = entry[keys == 'SubpixelPositions']
            subpix = entry['SubpixelPositions']
            if subpix == '0':
                subpix = [[1]]
            if subpix == '4-Point':
                subpix = [[4]]
            if subpix == '9-Point':
                subpix = [[9]]
            #in WFSS, SubpixelPositions will be either '4-Point' or '9-Point'
            #primary = entry[keys == 'PrimaryDithers']
            primary = entry['PrimaryDithers']
            if primary == '0':
                primary = [1]
            reps = np.int(subpix[0][0]) * np.int(primary[0])
            for key in keys:
                for j in range(reps):
                    expanded[key].append(indict[key][i])
        return expanded

    
    def base36encode(self,integer):
        chars, encoded = '0123456789abcdefghijklmnopqrstuvwxyz', ''

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            encoded = chars[remainder] + encoded

        return encoded.zfill(2)


    def get_pointing_info(self,file,propid):
        #read in information from APT's pointing file
        tar = []
        tile = []
        exp = []
        dith = []
        aperture = []
        targ1 = []
        targ2 = []
        ra = []
        dec = []
        basex = []
        basey = []
        dithx = []
        dithy = []
        v2 = []
        v3 = []
        idlx = []
        idly = []
        level = []
        type_str = []
        expar = []
        dkpar = []
        ddist = []
        observation_number = []
        visit_number = []
        visit_id = []
        visit_grp = []
        activity_id = []
        observation_label = []
        observation_id = []
        seq_id = []

        act_counter = 1
        with open(file) as f:
            for line in f:
                if len(line) > 1:
                    elements = line.split()
                    
                    #look for lines that give visit/observation numbers
                    if line[0:2] == '* ':
                        paren = line.rfind('(')
                        if paren == -1:
                            obslabel = line[2:]
                            obslabel = obslabel.strip()
                        else:
                            obslabel = line[2:paren-1]
                            obslabel = obslabel.strip()
                    if line[0:2] == '**':
                        v = elements[2]
                        obsnum,visitnum = v.split(':')
                        obsnum = str(obsnum).zfill(3)
                        visitnum = str(visitnum).zfill(3)
                        
                    try:
                        #skip the line at the beginning of each
                        #visit that gives info on the target,
                        #but is not actually an observation
                        #These lines have 'Exp' values of 0,
                        #while observations have a value of 1
                        #(that I've seen so far)
                        #
                        #Also, skip non-NIRCam lines. Check for NRC in aperture name
                        if ((np.int(elements[1]) > 0) & ('NRC' in elements[4])):
                            act = self.base36encode(act_counter)
                            activity_id.append(act)
                            observation_label.append(obslabel)
                            observation_number.append(obsnum)
                            visit_number.append(visitnum)
                            vid = str(propid)+visitnum+obsnum
                            visit_id.append(vid)
                            vgrp = '01'
                            visit_grp.append(vgrp)
                            seq = '1'
                            seq_id.append(seq)
                            tar.append(np.int(elements[0]))
                            tile.append(np.int(elements[1]))
                            exnum = str(elements[2]).zfill(5)
                            exp.append(exnum)
                            dith.append(np.int(elements[3]))

                            ap = elements[4]
                            if ('GRISMR_WFSS' in elements[4]):
                                ap = ap.replace('GRISMR_WFSS','FULL')
                            elif ('GRISMC_WFSS' in elements[4]):
                                ap = ap.replace('GRISMC_WFSS','FULL')
                            
                            aperture.append(ap)
                            targ1.append(np.int(elements[5]))
                            targ2.append(elements[6])
                            ra.append(elements[7])
                            dec.append(elements[8])
                            basex.append(elements[9])
                            basey.append(elements[10])
                            dithx.append(np.float(elements[11]))
                            dithy.append(np.float(elements[12]))
                            v2.append(np.float(elements[13]))
                            v3.append(np.float(elements[14]))
                            idlx.append(np.float(elements[15]))
                            idly.append(np.float(elements[16]))
                            level.append(elements[17])
                            type_str.append(elements[18])
                            expar.append(np.int(elements[19]))
                            dkpar.append(np.int(elements[20]))
                            ddist.append(np.float(elements[21]))
                            observation_id.append('V'+vid+'P00000000'+vgrp+seq+act)
                            act_counter += 1
                    except:
                        pass
                    
        pointing = {'exposure':exp, 'dither':dith, 'aperture':aperture,
                    'targ1':targ1, 'targ2':targ2, 'ra':ra, 'dec':dec,
                    'basex':basex, 'basey':basey, 'dithx':dithx,
                    'dithy':dithy, 'v2':v2, 'v3':v3, 'idlx':idlx,
                    'idly':idly, 'obs_label':observation_label,
                    'obs_num':observation_number,'visit_num':visit_number,
                    'act_id':activity_id,'visit_id':visit_id,'visit_group':visit_grp,
                    'sequence_id':seq_id,'observation_id':observation_id}
        return pointing


    def combine_dicts(self,dict1,dict2):
        #Now combine the dictionaries from the xml file and the pointing file
        combined = dict1.copy()
        combined.update(dict2)
        return combined


    def expand_for_detectors(self,dict):
        #Expand dictionary to have one line per detector, rather than the
        #one line per module that is in the input
        finaltab = {}
        for key in dict:
            finaltab[key] = []
        finaltab['detector'] = []
        for i in range(len(dict['PrimaryDithers'])):
            module = dict['Module'][i]
            if module == 'ALL':
                detectors = ['A1','A2','A3','A4','A5','B1','B2','B3','B4','B5']
            elif module == 'A':
                detectors = ['A1','A2','A3','A4','A5']
            elif module == 'B':
                detectors = ['B1','B2','B3','B4','B5']

            for key in dict:
                finaltab[key].extend(([dict[key][i]]*len(detectors)))
            finaltab['detector'].extend(detectors)
        return finaltab


    def ra_dec_update(self):
        #given the v2,v3 values in each entry, calculate RA,Dec

        #read in siaf
        distortionTable = ascii.read(self.siaf,header_start=1)
            
        aperture_ra = []
        aperture_dec = []
        for i in range(len(self.exposure_tab['Module'])):

            #first find detector
            #need ra,dec and v2,v3 pairs from entry
            #to calculate ra,dec at each detector's reference location
            detector = 'NRC' + self.exposure_tab['detector'][i]
            sub = self.exposure_tab['Subarray'][i]
            aperture = detector + '_' + sub
            pointing_ra = np.float(self.exposure_tab['ra'][i])
            pointing_dec = np.float(self.exposure_tab['dec'][i])
            pointing_v2 = np.float(self.exposure_tab['v2'][i])
            pointing_v3 = np.float(self.exposure_tab['v3'][i])
            pav3 = np.float(self.exposure_tab['pav3'][i])

            #calculate local roll angle
            local_roll = set_telescope_pointing.compute_local_roll(pav3,pointing_ra,
                                                                   pointing_dec,
                                                                   pointing_v2,
                                                                   pointing_v3)
            #create attitude_matrix    
            attitude_matrix = rotations.attitude(pointing_v2,pointing_v3,
                                                 pointing_ra,pointing_dec,local_roll)
            
            #find v2,v3 of the reference location for the detector
            aperture_v2,aperture_v3 = self.ref_location(distortionTable,aperture)
            
            #calculate RA, Dec of reference location for the detector
            ra,dec = rotations.pointing(attitude_matrix,aperture_v2,aperture_v3)
            #if dec < 0:
            #    print("ra, dec: {}, {}".format(ra,dec))
            #    print('attitude matrix {}'.format(attitude_matrix))
            #    print('aperture v2,v3: {}, {}'.format(aperture_v2.data,aperture_v3.data))
            #    print(pointing_v2,pointing_v3,pointing_ra,pointing_dec,local_roll,pav3)
            #    stop

            aperture_ra.append(ra)
            aperture_dec.append(dec)
        self.exposure_tab['ra_ref'] = aperture_ra
        self.exposure_tab['dec_ref'] = aperture_dec

        
    def ref_location(self,siaf,det):
        #find v2,v3 of detector reference location
        match = siaf['AperName'] == det
        if np.any(match) == False:
            print("Aperture name {} not found in input CSV file.".
                  format(aperture))
            sys.exit()
        v2 = siaf[match]['V2Ref']
        v3 = siaf[match]['V3Ref']
        return v2,v3

        
    def create_input_table(self):
        # Expand paths to full paths
        self.input_xml = os.path.abspath(self.input_xml)
        self.pointing_file = os.path.abspath(self.pointing_file)
        self.siaf = os.path.abspath(self.siaf)
        if self.output_csv is not None:
            self.output_csv = os.path.abspath(self.output_csv)
        if self.observation_table is not None:
            self.observation_table = os.path.abspath(self.observation_table)
            
        # Read in xml file
        tab = self.read_xml(self.input_xml)

        ascii.write(Table(tab), 'as_read_in.csv', format='csv', overwrite=True)
        
        # Expand the dictionary for multiple dithers. Expand such that there
        # is one entry in each list for each exposure, rather than one entry
        # for each set of dithers
        xmltab = self.expand_for_dithers(tab)

        #ascii.write(Table(xmltab), 'expand_for_dithers.csv', format='csv', overwrite=True)
        
        #read in the pointing file and produce dictionary
        pointing_tab = self.get_pointing_info(self.pointing_file,xmltab['ProposalID'][0])
        
        #combine the dictionaries
        obstab = self.combine_dicts(xmltab,pointing_tab)

        #ascii.write(Table(obstab), 'add_pointing_info.csv', format='csv', overwrite=True)

        #add epoch information
        #obstab = self.add_epochs(obstab)

        # add epoch and catalog information
        obstab = self.add_observation_info(obstab)
        
        # Expand for detectors. Create one entry in each list for each
        # detector, rather than a single entry for 'ALL' or 'BSALL'
        self.exposure_tab = self.expand_for_detectors(obstab)

        #ascii.write(Table(self.exposure_tab), 'expand_for_detectors.csv', format='csv', overwrite=True)
        
        # Calculate the correct V2,V3 and RA,Dec for each exposure/detector
        self.ra_dec_update()

        # Output to a csv file. 
        if self.output_csv is None:
            indir,infile = os.path.split(self.input_xml)
            self.output_csv = os.path.join(indir,'Observation_table_for_'+infile+'.csv')
        ascii.write(Table(self.exposure_tab), self.output_csv, format='csv', overwrite=True)
        print('Final csv exposure list written to {}'.format(self.output_csv))

        
    def read_observation_table(self,file):
        with open(file,'r') as infile:
            self.obstab = yaml.load(infile)
        

    def add_observation_info(self,intab):
        # Add information about each observation. Catalog names,
        # dates, PAV3 values, etc.
        self.read_observation_table(self.observation_table)

        onames = []
        for key1 in self.obstab:
            onames.append(self.obstab[key1]['Name'])
        onames = np.array(onames)
            
        obs_start = []
        obs_pav3 = []
        obs_sw_ptsrc = []
        obs_sw_galcat = []
        obs_sw_ext = []
        obs_sw_extscl = []
        obs_sw_extcent = []
        obs_sw_movptsrc = []
        obs_sw_movgal = []
        obs_sw_movext = []
        obs_sw_movconv = []
        obs_sw_solarsys = []
        obs_sw_bkgd = []
        obs_lw_ptsrc = []
        obs_lw_galcat = []
        obs_lw_ext = []
        obs_lw_extscl = []
        obs_lw_extcent = []
        obs_lw_movptsrc = []
        obs_lw_movgal = []
        obs_lw_movext = []
        obs_lw_movconv = []
        obs_lw_solarsys = []
        obs_lw_bkgd = []

        # This will allow the filter values
        # in the observation table to override
        # what comes from APT. This is useful for
        # WFSS observations, where you want to create
        # 'direct images' in various filters to hand
        # to the disperser software in order to create
        # simulated dispersed data. In that case, you
        # would need to create 'direct images' using
        # several filters inside the filter listed in APT
        # in order to get broadband colors to disperse.
        # If filter names are present in the observation
        # yaml file, then they will be used. If not, the
        # filters from APT will be kept
        obs_sw_filt = []
        obs_lw_filt = []

        for obs in intab['obs_label']:
            match = np.where(obs == onames)[0]
            if len(match) == 0:
                print("No valid epoch line found for observation {}".format(obs))
                print(type(obs))
                print(onames)
                sys.exit()
            else:
                #print('Matching {} from xml with {} from observation listfile'.format(obs,onames[match[0]]))
                obslist = self.obstab['Observation{}'.format(match[0]+1)]
                obs_start.append(obslist['Date'].strftime('%Y-%m-%d'))
                obs_pav3.append(obslist['PAV3'])
                obs_sw_ptsrc.append(obslist['SW']['PointSourceCatalog'])
                obs_sw_galcat.append(obslist['SW']['GalaxyCatalog'])
                obs_sw_ext.append(obslist['SW']['ExtendedCatalog'])
                obs_sw_extscl.append(obslist['SW']['ExtendedScale'])
                obs_sw_extcent.append(obslist['SW']['ExtendedCenter'])
                obs_sw_movptsrc.append(obslist['SW']['MovingTargetList'])
                obs_sw_movgal.append(obslist['SW']['MovingTargetSersic'])
                obs_sw_movext.append(obslist['SW']['MovingTargetExtended'])
                obs_sw_movconv.append(obslist['SW']['MovingTargetConvolveExtended'])
                obs_sw_solarsys.append(obslist['SW']['MovingTargetToTrack'])
                obs_sw_bkgd.append(obslist['SW']['BackgroundRate'])
                obs_lw_ptsrc.append(obslist['LW']['PointSourceCatalog'])
                obs_lw_galcat.append(obslist['LW']['GalaxyCatalog'])
                obs_lw_ext.append(obslist['LW']['ExtendedCatalog'])
                obs_lw_extscl.append(obslist['LW']['ExtendedScale'])
                obs_lw_extcent.append(obslist['LW']['ExtendedCenter'])
                obs_lw_movptsrc.append(obslist['LW']['MovingTargetList'])
                obs_lw_movgal.append(obslist['LW']['MovingTargetSersic'])
                obs_lw_movext.append(obslist['LW']['MovingTargetExtended'])
                obs_lw_movconv.append(obslist['LW']['MovingTargetConvolveExtended'])
                obs_lw_solarsys.append(obslist['LW']['MovingTargetToTrack'])
                obs_lw_bkgd.append(obslist['LW']['BackgroundRate'])

                # Override filters if given
                try:
                    obs_sw_filt.append(obslist['SW']['Filter'])
                except:
                    pass
                try:
                    obs_lw_filt.append(obslist['LW']['Filter'])
                except:
                    pass

        intab['epoch_start_date'] = obs_start
        intab['pav3'] = obs_pav3
        intab['sw_ptsrc'] = obs_sw_ptsrc
        intab['sw_galcat'] = obs_sw_galcat
        intab['sw_ext'] = obs_sw_ext
        intab['sw_extscl'] = obs_sw_extscl
        intab['sw_extcent'] = obs_sw_extcent
        intab['sw_movptsrc'] = obs_sw_movptsrc
        intab['sw_movgal'] = obs_sw_movgal
        intab['sw_movext'] = obs_sw_movext
        intab['sw_movconv'] = obs_sw_movconv
        intab['sw_solarsys'] = obs_sw_solarsys
        intab['sw_bkgd'] = obs_sw_bkgd
        intab['lw_ptsrc'] = obs_lw_ptsrc
        intab['lw_galcat'] = obs_lw_galcat
        intab['lw_ext'] = obs_lw_ext
        intab['lw_extscl'] = obs_lw_extscl
        intab['lw_extcent'] = obs_lw_extcent
        intab['lw_movptsrc'] = obs_lw_movptsrc
        intab['lw_movgal'] = obs_lw_movgal
        intab['lw_movext'] = obs_lw_movext
        intab['lw_movconv'] = obs_lw_movconv
        intab['lw_solarsys'] = obs_lw_solarsys
        intab['lw_bkgd'] = obs_lw_bkgd

        # Here we override the filters read from APT
        # if they are given in the observation yaml file
        if len(obs_sw_filt) > 0:
            intab['ShortFilter'] = obs_sw_filt
        if len(obs_lw_filt) > 0:
            intab['LongFilter'] = obs_lw_filt
        return intab
    
        
    def add_epochs(self,intab):
        #add information on the epoch of each observation
        #if the user entered a list of epochs, read that in
        default_date = '2020-10-14'
        
        if self.epoch_list is not None:
            epochs = ascii.read(self.epoch_list,header_start=0,data_start=1)
        else:
            epochs = Table()
            epochs['observation'] = intab['obs_label']
            epochs['date'] = ['2018-10-14'] * len(intab['obs_label'])
            epochs['pav3'] = [0.] * len(intab['obs_label'])
            
        #insert epoch info for each observation into dictionary
        epoch_start = []
        epoch_pav3 = []
        for obs in intab['obs_label']:
            match = obs == epochs['observation'].data
            if np.sum(match) == 0:
                print("No valid epoch line found for observation {}".format(obs))
                print(epochs['observation'].data)
                epoch_start.append(default_date)
                epoch_pav3.append(0.)
            else:
                epoch_start.append(epochs['date'][match].data[0])
                epoch_pav3.append(epochs['pav3'][match].data[0])
        intab['epoch_start_date'] = epoch_start
        intab['pav3'] = epoch_pav3
        return intab
        

    def dict_lengths(self,dict):
        minlength = 99999999
        maxlength = 0
        for key in dict:
            ll = len(dict[key])
            if ll > maxlength:
                maxlength = ll
            if ll < minlength:
                minlength = ll
        return minlength,maxlength
    

    def add_options(self,parser=None,usage=None):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,description='Simulate JWST ramp')
        parser.add_argument("input_xml",help='XML file from APT describing the observations.')
        parser.add_argument("pointing_file",help='Pointing file from APT describing observations.')
        parser.add_argument("siaf",help='Name of CSV version of SIAF')
        parser.add_argument("--output_csv",help="Name of output CSV file containing list of observations.",default=None)
        parser.add_argument("--observation_table",help='Ascii file containing a list of observations, times, and roll angles, catalogs',default=None)
        return parser
    

if __name__ == '__main__':

    usagestring = 'USAGE: apt_inputs.py NIRCam_obs.xml NIRCam_obs.pointing SIAF_March2017.csv'

    input = AptInput()
    parser = input.add_options(usage = usagestring)
    args = parser.parse_args(namespace=input)
    input.create_input_table()
