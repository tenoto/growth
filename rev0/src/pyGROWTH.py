#!/usr/bin/env python

__pyname__  = 'pyGROWTH.py'
__version__ = 'rev0.01'
__date__    = '2015-11-06'
__project__ = 'GROWTH collaboration'
__author__  = 'Teruaki Enoto'
__email__   = 'teruaki.enoto@gmail.com'

import os 
import sys 
import yaml 
import ROOT 
#import array
import pyfits
import numpy 
import time 

DEBUG = True

TOOLBAR_WIDTH = 40
A4_WIDTH   = 842
#A4_WIDTH_2 = 800
A4_HEIGHT  = 595 

FPGA_CLOCK_PERIOD = 10e-9 # 10 ns, 100 MHz 
FPGA_CLOCK_BIT = 40 # bit 

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)

class GrowthFitsFile():
	def __init__(self, fitsfile, configyaml='./config/pyGROWTH_config.yaml'):
		sys.stdout.write("==========================================\n")
		sys.stdout.write("pyGROWTH: GrowthFitsFile initialization...\n")
		self.fitsfile = fitsfile

		if not os.path.exists(self.fitsfile):
			sys.stderr.write("Error: file, %s, does not exists.\n" % self.fitsfile)
			quit()
		self.fitsname = os.path.basename(self.fitsfile)
		self.basename = os.path.splitext(self.fitsname)[0]
		try:
			self.hdu  = pyfits.open(self.fitsfile)
		except:
			sys.stderr.write("Error: fits file, %s, can not be opended.\n" % self.fitsfile)
			quit()			
		sys.stdout.write("Input fitsfile %s is successfully loaded.\n" % self.fitsfile) 

		self.configyaml = configyaml
		if not os.path.exists(self.configyaml):
			sys.stderr.write("Error: file, %s, does not exists.\n" % self.configyaml)
			quit()			
		f = open(self.configyaml)
		self.config = yaml.load(f)
		f.close()
		sys.stdout.write("Config file %s is successfully loaded.\n" % self.configyaml) 		

	def setResultDictionary(self):
		sys.stdout.write('[setResultDictionary]\n')		

		self.result = {}
		history_flag = True
		for keyword in self.hdu['EVENTS'].header:
			if keyword == 'HISTORY':
				if history_flag:
					for line in self.hdu['EVENTS'].header[keyword]:
						cols = line.split('YAML-- ')[-1].split(':')
						if len(cols) == 2:
							try:
								self.result['DAQCONFIG_%s' % cols[0]] = int(cols[1])
							except:
								self.result['DAQCONFIG_%s' % cols[0]] = cols[1]
				history_flag = False
				continue
			self.result['HEADER_%s' % keyword] = self.hdu['EVENTS'].header[keyword]

		self.result['TOTAL_COUNTS'] = int(self.hdu['EVENTS'].header['NAXIS2'])
		self.result['RAW_RATE']     = float(self.hdu['EVENTS'].header['NAXIS2']) / float(self.result['HEADER_EXPOSURE'])

	def makeOutputDirectory(self, overwrite=True):
		sys.stdout.write('[makeOutputDirectory]\n')		

		self.outdir = '%s/%s' % (self.config['outroot'],self.basename)
		if overwrite:
			os.system('rm -rf %s' % self.outdir)
		os.system('mkdir -p %s' % self.outdir)

	def convertToROOT(self):
		sys.stdout.write('[convertToROOT] Converting to ROOT.\n')		

		self.rootfile = '%s/%s.root' % (self.outdir, self.basename)
		self.tfile = ROOT.TFile(self.rootfile, 'recreate')
		self.ttree = ROOT.TTree('EVENTS','EVENTS')
		timeTag = numpy.array([0],dtype='int64')
		# http://oxon.hatenablog.com/entry/20100506/1273102651
		self.ttree.Branch('timeTag',timeTag,'timeTag/l')

		timeTagElapse = numpy.array([0],dtype='int64')
		self.ttree.Branch('timeTagElapse',timeTagElapse,'timeTagElapse/l')		

		deltaTimeTag = numpy.array([0],dtype='int64')
		self.ttree.Branch('deltaTimeTag',deltaTimeTag,'deltaTimeTag/l')				

		deltaTriggerCount = numpy.array([0],dtype='int64')
		self.ttree.Branch('deltaTriggerCount',deltaTriggerCount,'deltaTriggerCount/l')					

		phaMax = numpy.array([0],dtype='int32')
		self.ttree.Branch('phaMax',phaMax,'phaMax/I')	

		baseline = numpy.array([0],dtype='int32')
		self.ttree.Branch('baseline',baseline,'baseline/I')			

		sys.stdout.write("[%s]" % (" " * TOOLBAR_WIDTH))
		sys.stdout.flush()
		sys.stdout.write("\b" * (TOOLBAR_WIDTH+1)) # return to start of line, after '['

		naxis2 = self.hdu['EVENTS'].header['NAXIS2']
		TOOLBAR_BASE = int(naxis2 / TOOLBAR_WIDTH)

		timeTag0 = self.hdu['EVENTS'].data['timeTag'][0]
		timeTagPrev = timeTag0

		triggerCount0 = self.hdu['EVENTS'].data['triggerCount'][0]
		triggerCountPrev = triggerCount0

		evtnum = 0
		for evt in self.hdu['EVENTS'].data:
			timeTag[0] = evt['timeTag']
			
			timeTagElapse[0] = evt['timeTag'] - timeTag0			
			if evt['timeTag'] < timeTagPrev:
				timeTagElapse[0] += 2**FPGA_CLOCK_BIT

			deltaTimeTag[0] =  evt['timeTag'] - timeTagPrev
			timeTagPrev = evt['timeTag']

			deltaTriggerCount[0] = evt['triggerCount'] - triggerCountPrev
			triggerCountPrev = evt['triggerCount']

			phaMax[0] = evt['phaMax']
			baseline[0] = evt['baseline']

			self.ttree.Fill()

		 	evtnum += 1

			if evtnum % TOOLBAR_BASE == 0:
				sys.stdout.write("-")
		 		sys.stdout.flush()

		sys.stdout.write("\n")
		sys.stdout.write("Conversion is completed!\n")		

	def DrawLightCurve(self):
		sys.stdout.write('[DrawLightCurve]\n')		

		name = '%s_rawlc' % self.basename
		nbin = int(int(self.result['HEADER_EXPOSURE'])/self.config['LC_TBIN'])

		can = ROOT.TCanvas('can','can',A4_WIDTH,A4_HEIGHT)
		can.SetMargin(0.1,0.03,0.13,0.07)

		self.ttree.Draw("timeTagElapse*%.1e>>%s(%d,0,%.1f)" % 
			(FPGA_CLOCK_PERIOD,name,nbin,self.result['HEADER_EXPOSURE']))
		self.th1_rawlc = self.tfile.Get(name)

		self.th1_rawlc.Fit('pol0','','',self.config['LC_FIT_XMIN'],self.result['HEADER_EXPOSURE'])
		self.th1_rawlc_fit = self.th1_rawlc.GetFunction('pol0')
		self.th1_rawlc_fit.SetLineColor(ROOT.kRed)

		self.th1_rawlc.SetTitle('%s;Time (sec) %.1f-s bin;Counts/bin' % (name,self.config['LC_TBIN']))
		self.th1_rawlc.SetMinimum(self.config['LC_PLOT_YMIN'])
		self.th1_rawlc.SetMaximum(self.config['LC_PLOT_YMAX'])
		self.th1_rawlc.Draw('hist e')
		self.th1_rawlc_fit.Draw('same')

		self.result['LC_FIT_NDF'] = self.th1_rawlc_fit.GetNDF()
		self.result['LC_FIT_CHI2'] = self.th1_rawlc_fit.GetChisquare()		
		self.result['LC_FIT_MEAN'] = self.th1_rawlc_fit.GetParameter(0)		
		self.result['LC_FIT_MEAN_ERROR'] = self.th1_rawlc_fit.GetParError(0)

		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)		

		name = '%s_rawlc_rate' % self.basename		
		#self.th1_rawlc_rate = self.th1_rawlc.Clone()
		self.th1_rawlc_rate = ROOT.TH1D(name,name,nbin,0.,self.result['HEADER_EXPOSURE'])
		#self.th1_rawlc_rate.SetName(name)
		for i in range(0,self.th1_rawlc.GetNbinsX()+1):
			self.th1_rawlc_rate.SetBinContent(i,self.th1_rawlc.GetBinContent(i)/self.config['LC_TBIN'])		
			self.th1_rawlc_rate.SetBinError(i,self.th1_rawlc.GetBinError(i)/self.config['LC_TBIN'])		
		self.th1_rawlc_rate.Fit('pol0')
		self.th1_rawlc_rate_fit = self.th1_rawlc_rate.GetFunction('pol0')
		self.th1_rawlc_rate_fit.SetLineColor(ROOT.kRed)
		self.result['LC_FIT_RATE'] = self.th1_rawlc_rate_fit.GetParameter(0)		
		self.result['LC_FIT_RATE_ERROR'] = self.th1_rawlc_rate_fit.GetParError(0)

		self.th1_rawlc_rate.SetTitle('%s;Time (sec) %.1f-s bin;Counts/sec' % (name,self.config['LC_TBIN']))		
		self.th1_rawlc_rate.SetMaximum(float(self.config['LC_PLOT_YMAX'])/float(self.config['LC_TBIN']))
		self.th1_rawlc_rate.Draw('hist e')
		self.th1_rawlc_rate_fit.Draw('same')
		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)	

		name = '%s_lchist' % self.basename
		self.th1_lchist = ROOT.TH1I(name,name,
			self.config['LCHIST_NBIN'],0.0,int(self.config['LC_PLOT_YMAX']))
		for i in range(0,self.th1_rawlc.GetNbinsX()+1):
			if self.th1_rawlc.GetBinContent(i) > 0:
				self.th1_lchist.Fill(self.th1_rawlc.GetBinContent(i))

		#print self.result['RAW_RATE']*self.config['LC_TBIN']
		ROOT.gROOT.GetFunction('gaus').SetParameter(1,self.result['RAW_RATE']*self.config['LC_TBIN'])
		#ROOT.gROOT.GetFunction('gaus').SetParLimits(1,0.1*self.result['RAW_RATE']*self.config['LC_TBIN'],0.9*self.result['RAW_RATE']*self.config['LC_TBIN'])
		self.th1_lchist.Fit('gaus',"","",0.1*self.result['RAW_RATE']*self.config['LC_TBIN'],1.9*self.result['RAW_RATE']*self.config['LC_TBIN'])
		#self.th1_lchist.Fit('gaus')
		self.th1_lchist_fit = self.th1_lchist.GetFunction('gaus')
		self.th1_lchist_fit.SetLineColor(ROOT.kRed)
		self.th1_lchist_fit.SetLineWidth(1)

		self.result['LCHIST_NDF']  = self.th1_lchist_fit.GetNDF()
		self.result['LCHIST_CHI2'] = self.th1_lchist_fit.GetChisquare()
		self.result['LCHIST_CONSTANT']       = self.th1_lchist_fit.GetParameter(0)
		self.result['LCHIST_CONSTANT_ERROR'] = self.th1_lchist_fit.GetParError(0)
		self.result['LCHIST_MEAN']       = self.th1_lchist_fit.GetParameter(1)
		self.result['LCHIST_MEAN_ERROR'] = self.th1_lchist_fit.GetParError(1)		
		self.result['LCHIST_SIGMA']       = self.th1_lchist_fit.GetParameter(2)
		self.result['LCHIST_SIGMA_ERROR'] = self.th1_lchist_fit.GetParError(2)				

		self.th1_lchist.SetTitle('%s;Counts/sec;Number' % name)
		self.th1_lchist.Draw('hist e')
		self.th1_lchist_fit.Draw('same')
		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)	

		can.Closed()

	def DrawDeltaTimeTag(self):
		sys.stdout.write('[DrawDeltaTimeTag]\n')		

		name = '%s_deltat' % self.basename
		nbin = int(float(self.config['DELTAT_PLOT_XMAX'])/float(self.config['DELTAT_TBIN']))

		can = ROOT.TCanvas('can','can',A4_WIDTH,A4_HEIGHT)
		can.SetMargin(0.1,0.03,0.13,0.07)
		can.SetLogy()	

		self.ttree.Draw("deltaTimeTag*%.1e>>%s(%d,0,%.1f)" % 
			(FPGA_CLOCK_PERIOD,name,nbin,self.config['DELTAT_PLOT_XMAX']))
		self.th1_deltat = self.tfile.Get(name)

		self.th1_deltat.Fit('expo')
		self.th1_deltat_fit = self.th1_deltat.GetFunction('expo')
		self.th1_deltat_fit.SetLineColor(ROOT.kRed)
		self.th1_deltat_fit.SetLineWidth(1)		

		self.th1_deltat.SetTitle('%s;Time (sec) %.1f-s bin;Counts/bin' % (name,self.config['DELTAT_TBIN']))
		self.th1_deltat.GetXaxis().SetLimits(0,self.config['DELTAT_PLOT_XMAX'])
		self.th1_deltat.GetXaxis().SetRangeUser(0,self.config['DELTAT_PLOT_XMAX'])

		self.th1_deltat.Draw('hist e')
		self.th1_deltat_fit.Draw('same')

		self.result['DELTAT_FIT_DOF']   = self.th1_deltat_fit.GetNDF()
		self.result['DELTAT_FIT_CHI2']  = self.th1_deltat_fit.GetChisquare()
		self.result['DELTAT_FIT_CONST'] = self.th1_deltat_fit.GetParameter(0)
		self.result['DELTAT_FIT_CONST_ERROR'] = self.th1_deltat_fit.GetParError(0)
		self.result['DELTAT_FIT_SLOPE'] = self.th1_deltat_fit.GetParameter(1)
		self.result['DELTAT_FIT_SLOPE_ERROR'] = self.th1_deltat_fit.GetParError(1)		

		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)		

		can.Closed()

	def DrawDeltaTriggerCount(self):
		sys.stdout.write('[DrawDeltaTriggerCount]\n')		
		name = '%s_deltaTrigCnt' % self.basename
		nbin = int(self.config['DELTATRIGGERCOUNT_PLOT_XMAX'])

		can = ROOT.TCanvas('can','can',A4_WIDTH,A4_HEIGHT)
		can.SetMargin(0.1,0.03,0.13,0.07)
		can.SetLogy()	

		self.ttree.Draw("deltaTriggerCount>>%s(%d,%.1f,%.1f)" % (
			name, nbin, 
			self.config['DELTATRIGGERCOUNT_PLOT_XMIN'],
			self.config['DELTATRIGGERCOUNT_PLOT_XMAX']))
		self.th1_deltaTriggerCount = self.tfile.Get(name)

		self.th1_deltaTriggerCount.Draw('hist e')
		#self.th1_deltaTriggerCount.Draw('e same')

		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)		

		can.Closed()

	def DrawSpectrum(self):
		sys.stdout.write('[DrawSpectrum]\n')

		name = '%s_spec' % self.basename
		nbin = int(self.config['SPEC_PLOT_XMAX']-self.config['SPEC_PLOT_XMIN'])

		can = ROOT.TCanvas('can','can',A4_WIDTH,A4_HEIGHT)
		can.SetMargin(0.1,0.03,0.13,0.07)
		can.SetLogy()	

		self.ttree.Draw("phaMax-baseline>>%s(%d,%.1f,%.1f)" % (
			name, nbin, 
			self.config['SPEC_PLOT_XMIN'],
			self.config['SPEC_PLOT_XMAX']))
		self.th1_spec = self.tfile.Get(name)
		self.th1_spec.SetTitle('%s;phaMax-baseline (ch);Counts/bin' % name)

		self.th1_spec.Draw('hist e')

		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)		

		name = '%s_normspec' % self.basename		
		#self.th1_normalized_spec = self.th1_spec.Clone()
		self.th1_normalized_spec = ROOT.TH1D(name,name,nbin,self.config['SPEC_PLOT_XMIN'],self.config['SPEC_PLOT_XMAX'])
		for i in range(0,self.th1_spec.GetNbinsX()+1):
			self.th1_normalized_spec.SetBinContent(i,self.th1_spec.GetBinContent(i)/self.result['HEADER_EXPOSURE'])		
			self.th1_normalized_spec.SetBinError(i,self.th1_spec.GetBinError(i)/self.result['HEADER_EXPOSURE'])		

		self.th1_normalized_spec.SetTitle('%s;phaMax-baseline (ch);Counts/sec' % name)
		self.th1_normalized_spec.Draw('hist e')
		outpdf = '%s/%s.pdf' % (self.outdir, name)
		can.Print(outpdf)	

		can.Closed()

	def close(self):
		sys.stdout.write('[close]\n')		

		self.outyaml = '%s/%s.yaml' % (self.outdir, self.basename)
		f = open(self.outyaml,'w')
		f.write(yaml.dump(self.result,
			default_flow_style=False, default_style=''))
		f.close()	

		try:
			print self.tfile.Write()
			if self.hdu['EVENTS'].header['NAXIS2'] > 1e+6:
				time.sleep(5)
			elif self.hdu['EVENTS'].header['NAXIS2'] > 1e+5:
				time.sleep(4)
			elif self.hdu['EVENTS'].header['NAXIS2'] > 1e+4:
				time.sleep(1)				
			else:
				time.sleep(1)	

			self.tfile.Close("R")
			sys.stdout.write('File is written to %s\n' % self.rootfile)
			sys.stdout.write('Finished!\n')
			return 0
		except:
			sys.stderr.write('Writing file was failed. Please try again.\n')					
			return 1 

	def run(self):
		self.makeOutputDirectory()
		self.setResultDictionary()
		self.convertToROOT()
		self.setResultDictionary()
		self.DrawLightCurve()
		self.DrawDeltaTimeTag()
		self.DrawDeltaTriggerCount()
		self.DrawSpectrum()
		self.close()

if len(sys.argv) != 2:
	sys.stderr.write('%s fitsfile \n' % sys.argv[0])
	quit()		

fits = GrowthFitsFile(sys.argv[1])		
fits.run()

