import os
import subprocess
import ROOT

ROOT.gErrorIgnoreLevel=ROOT.kWarning # to remove printing pdf page message

allfiles=os.listdir('./')

# model='ILD_l5_o1_v02_'
model='ILD_s5_o1_v02_'
setuplab='tauFind_ALLEVT_'+model
# setuplab='tauFind_'+model

prelab='results_'+setuplab

# pol='unpol'
pol='eL80pR30'
# pol='eR80pL30'
# pol='eL80pL30'
# pol='eR80pR30'


print '**********************************************'
print model
print '**********************************************'

logs=[]
for ff in allfiles:
    if prelab in ff:
        gg=ff.split('_')
        proc=gg[len(gg)-2]
        if proc.startswith('I'):
            logs.append(ff)

logs.sort()

#--------------------------

saveFitHistos=True
plotHistos=True

if plotHistos or saveFitHistos:
    figdir='./figs/'+model+pol+'/'
    if not os.path.isdir(figdir):
        os.makedirs(figdir)

fithistos=0
if saveFitHistos:
    fithistos=ROOT.TFile(figdir+"fithistos.root",'recreate')

#--------------------------

totlumi=4000. # in fb-1

elePol=0.
posPol=0.
lumi=0.
if pol=='eL80pR30':
    elePol=-0.8
    posPol=0.3
    lumi=totlumi*0.4
elif pol=='eR80pL30':
    elePol=0.8
    posPol=-0.3
    lumi=totlumi*0.4
elif pol=='eL80pL30':
    elePol=-0.8
    posPol=-0.3
    lumi=totlumi*0.1
elif pol=='eR80pR30':
    elePol=0.8
    posPol=0.3
    lumi=totlumi*0.1
else:
    pol='unpol'
    lumi=totlumi

ttf=ROOT.TText()
ttf.SetTextFont(62)
ttf.SetTextSize(0.06)
labPosX=0.12
labPosY=0.92
ILDlab="ILD preliminary"

ttf2=ROOT.TText()
ttf2.SetTextSize(0.03)
ttf2.SetTextAngle(-90.)
lab2PosX=0.92
lab2PosY=0.6
beamLab = model+pol+'_L'+str(int(lumi))



# get the polarisation weights
wt_eL = (1 + abs(elePol))/2
wt_eR = 1. - (1 + abs(elePol))/2
if elePol>0:
    temp=wt_eL
    wt_eL=wt_eR
    wt_eR=temp
wt_pL =  (1 + abs(posPol))/2
wt_pR =  1. - (1 + abs(posPol))/2
if posPol>0:
    temp=wt_pL
    wt_pL=wt_pR
    wt_pR=temp

print '**********************************************'
print '    elePol=', elePol, 'posPol=', posPol, 'lumi=', lumi
print '**********************************************'

print '%50s' % 'process',
print '%10s' % 'XSEC',
print '%10s' % 'nMCorig',
print '%10s' % 'nMCsel',
print '%10s' % 'nExpect(k)'


nclass=7

totalSignal=0
totalBackground_tt=0
totalBackground_2f=0
totalBackground_4f=0
totalBackground_oth=0

splitSignal={}
splitBackground_tt={}
splitBackground_2f={}
splitBackground_4f={}
splitBackground_oth={}


# sigweight={}
procweight={}

cutnames=['nMC']
cutevents={}

ngroupOrigProc={}
#ngroupOrigProc['signal']={}
#ngroupOrigProc['bg2f']={}
#ngroupOrigProc['bg4f']={}
#ngroupOrigProc['bgoth']={}
procgroupnames=['signal', 'bg2f', 'bg4f'] # , 'bgoth']

subclassnames=[]

for logfilename in logs:

    signal = '2f_Z_leptonic' in logfilename

    if signal:
        procgroup=procgroupnames[0]
    elif '2f' in logfilename:
        procgroup=procgroupnames[1]
    elif '4f' in logfilename:
        procgroup=procgroupnames[2]
    else:
        print "ERROR unable to assign to process group!", logfilename
        break

    pname=logfilename.replace(prelab,'').replace('.txt','')

    gg=logfilename.split('_')
    processID=gg[len(gg)-2]
    
    ret=subprocess.check_output(['grep', processID, 'input_data_info.csv'] )
    xsec = float( (ret.split(';')[15]).strip('"') )

    polE = gg[len(gg)-1].split('.')[0]
    polP = gg[len(gg)-1].split('.')[1]

    polWt=1
    if polE=='eL':
        polWt=polWt*wt_eL
    elif polE=='eR':
        polWt=polWt*wt_eR
    else :
        print 'unknown electron polarisation', polE
        continue
        
    if polP=='pL':
        polWt=polWt*wt_pL
    elif polP=='pR':
        polWt=polWt*wt_pR
    else :
        print 'unknown positron polarisation', polP
        continue

    polxsec = xsec*polWt

    logfile=open(logfilename,'read')
    lines=logfile.readlines()
    insummary=False

    nsplitsplitsel={} # separated over sample subclasses

    norig=0
    for line in lines:
        ff=line.split()

        if len(subclassnames)==0 and 'TOTAL' in line:
            for i in range(nclass):
                subclassnames.append(ff[i])

        lab=ff[0].strip()
        if lab=='_nOrig' or lab=='_nPresel' or lab=='_nTwoSeeds' or '_nCumulSel' in lab:
            if lab not in cutnames:
                cutnames.append(lab)
            nsplitsplitsel[lab]=[]

            # calculate total evts in sample, and sample weight
            if lab=='_nOrig':
                nsplitsplitsel['nMC']=[]
                for i in range(nclass):
                    nsplitsplitsel['nMC'].append( int(ff[i+1]) )
                    norig=norig+int(ff[i+1])
                procweight[pname]=lumi*polxsec/norig

            for i in range(nclass):
                nsplitsplitsel[lab].append( int( procweight[pname]*int(ff[i+1]) ) ) # the number of exp events
                

    #print cutnames
    #for cut in cutnames:
    #    print cut , nsplitsplitsel[cut]

    if procgroup not in ngroupOrigProc.keys():
        ngroupOrigProc[procgroup] = nsplitsplitsel
    else:
        for cut in nsplitsplitsel.keys():
            for i in range(len( nsplitsplitsel[cut] ) ):
                ngroupOrigProc[procgroup][cut][i]=ngroupOrigProc[procgroup][cut][i]+nsplitsplitsel[cut][i]
    
# print ngroupOrigProc.keys()
# for proc in ngroupOrigProc.keys():
#     print proc, ngroupOrigProc[proc]


print '\n \n'
print 'COMPLETE TABLE'
print '%20s' % (''),
for pgn in procgroupnames:
    print '%25s' % '',
    print '%24s' % pgn,
    print '%25s' % '',
    print ' : ',
print ''

print '%20s' % (''),
for i in range(len(procgroupnames)):
    for subn in subclassnames:
        print '%10s' % subn,
    print ' : ',
print ''

for cut in cutnames:
    print '%20s' % (cut),
    for proc in procgroupnames:
        for i in range(nclass):
            print '%10d' % ngroupOrigProc[proc][cut][i],
        print ' : ',
    print ''


print '\n \n'
print 'REDUCED TABLE'
print '%20s' % (''),
print '%33s' % procgroupnames[0],
print '%20s' % (''),
print ' : ',
print '%24s' % procgroupnames[1],
print '%18s' % (''),
print ' : ',
print '%24s' % procgroupnames[2]
print ''
print '%20s' % (''),
print '%10s' % subclassnames[0],
print '%10s' % subclassnames[1],
print '%10s' % subclassnames[2],
print '%10s' % subclassnames[3],
print '%10s' % 'OTHER',
for i in range(2):
    print ' : ',
    print '%10s' % '0-TAU',
    print '%10s' % '1-TAU',
    print '%10s' % '2-TAU',
    print '%10s' % '>2-TAU',
print ''
print ''

for cut in cutnames:
    print '%20s' % (cut),
    for proc in procgroupnames:
        thissplit=[]
        if proc=='signal':
            thissplit.append( ngroupOrigProc[proc][cut][0] )
            thissplit.append( ngroupOrigProc[proc][cut][1] )
            thissplit.append( ngroupOrigProc[proc][cut][2] )
            thissplit.append( ngroupOrigProc[proc][cut][3] )
            thissplit.append( ngroupOrigProc[proc][cut][4]+ngroupOrigProc[proc][cut][5]+ngroupOrigProc[proc][cut][6] )
        else:
            thissplit.append( ngroupOrigProc[proc][cut][4] )
            thissplit.append( ngroupOrigProc[proc][cut][5] )
            thissplit.append( ngroupOrigProc[proc][cut][0]+ngroupOrigProc[proc][cut][1]+ngroupOrigProc[proc][cut][2]+ngroupOrigProc[proc][cut][3] )
            thissplit.append( ngroupOrigProc[proc][cut][6] )

        for gt in thissplit:
            print '%10d' % gt ,
        print ' : ',
    print ''



print '\n \n'
print 'REDUCED TABLE2'
print '%20s' % (''),
print '%20s' % procgroupnames[0],
print '%11s' % (''),
print ' : ',
print '%10s' % procgroupnames[1],
print ' :',
print '%10s' % procgroupnames[2]
print ''
print '%20s' % (''),
print '%10s' % ('eff'),
print '%10s' % ('nevt'),
print '%10s' % ('OTH'),
print ''

for cut in cutnames:
    if cut=='nMC':
        continue

    print '%20s' % (cut),
    for proc in procgroupnames:
        thissplit=[]
        if proc=='signal':
            thissplit.append( float(ngroupOrigProc[proc][cut][0])/float(ngroupOrigProc[proc]['_nOrig'][0]) )
            thissplit.append( ngroupOrigProc[proc][cut][0] )
            thissplit.append( ngroupOrigProc[proc][cut][1]+ngroupOrigProc[proc][cut][2]+ngroupOrigProc[proc][cut][3]+ngroupOrigProc[proc][cut][4]+ngroupOrigProc[proc][cut][5]+ngroupOrigProc[proc][cut][6] )

            print '%10.1f' % (100.*thissplit[0]),
            for i in range(1, len(thissplit) ):
                print '%10.1f' % float(thissplit[i]/1000.),

        else:
            thissplit.append( ngroupOrigProc[proc][cut][0]+ngroupOrigProc[proc][cut][1]+ngroupOrigProc[proc][cut][2]+ngroupOrigProc[proc][cut][3]+ ngroupOrigProc[proc][cut][4]+ ngroupOrigProc[proc][cut][5]+ ngroupOrigProc[proc][cut][6] )

            print '%10.1f' % float(thissplit[0]/1000.),
        print ' : ',
    print ''



if plotHistos or saveFitHistos:

    # group into classes for plotting
    classes={}
    # cols=[1,2,4,6,8,9,13,14,15]
    cols=[ROOT.kRed, ROOT.kOrange, ROOT.kMagenta+1, ROOT.kBlue-3, ROOT.kAzure+7, ROOT.kCyan+3, ROOT.kCyan-9]
    classCols={}
    icol=1
    for logfilename in logs:
        basename=logfilename.replace(prelab,'').replace('.txt','')
        classname=basename[0:len(basename)-14]
        if classname not in classes:
            classes[classname]=[]
            classCols[classname]=cols[icol%len(cols)]
            icol=icol+1
        classes[classname].append(basename)

    # print classCols

    histonames=[]
    rootfiles={}

    for cl in classes.keys():
        rootfiles[cl]={}
        for subcl in classes[cl]:
            # get list of histograms we want to plot
            rootname=setuplab+subcl+'.root'

            print rootname

            rootfile=ROOT.TFile(rootname, 'read')
            if rootfile.GetNkeys() == 0:
                rootfile.Close()
            else:
                rootfiles[cl][subcl]=rootfile
                if len(histonames)==0:
                    allkeys=rootfile.GetListOfKeys()
                    for key in allkeys:
                        obj=key.ReadObj()
                        if obj.ClassName()=='TH1F' or obj.ClassName()=='TH1D':
                            histonames.append(obj.GetName())

    orderedClasses=[]
    if '2f_Z_leptonic' in classes.keys():
        orderedClasses.append('2f_Z_leptonic')
    dddd=classes.keys()
    dddd.sort()
    for cl in  dddd:
        if cl != '2f_Z_leptonic':
            orderedClasses.append(cl)

    histonames2=[]
    for histoname in histonames:
        if ( plotHistos and 'sample0' in histoname ) or \
                ( saveFitHistos and 'sample0_SEL_rec_' in histoname and '_pol' in histoname ):
            histonames2.append(histoname)
            
            # ) or ( doFits and 'sample0' in histoname and 'SEL_rec_rho_pol' in histoname):

    grouping={}
    grouping['2f_Z_leptonic']={}
    grouping['2f_Z_leptonic'][0]=0  # 2f 2tau sig
    grouping['2f_Z_leptonic'][1]=1  # 2f 2tau bkg
    grouping['2f_Z_leptonic'][2]=1
    grouping['2f_Z_leptonic'][3]=1
    grouping['2f_Z_leptonic'][4]=2  # 2f not 2tau 
    grouping['2f_Z_leptonic'][5]=2
    grouping['2f_Z_leptonic'][6]=2

    grouping['2f_Z_hadronic']={}
    grouping['2f_Z_bhabhag']={}
    for i in range(nclass):
        grouping['2f_Z_bhabhag'][i]=2
        grouping['2f_Z_hadronic'][i]=2

    grouping['4f']={}
    grouping['4f'][0]=5  # 4f 2T
    grouping['4f'][1]=5
    grouping['4f'][2]=5
    grouping['4f'][3]=5
    grouping['4f'][4]=3  # 4f 0T
    grouping['4f'][5]=4  # 4f 1T
    grouping['4f'][6]=6  # 4f >2T

    grouplabels={}
    grouplabels[0]="2-tau signal"
    grouplabels[1]="2-tau bkg"
    grouplabels[2]="other 2-f"
    grouplabels[3]="4-f (0 tau)"
    grouplabels[4]="4-f (1 tau)"
    grouplabels[5]="4-f (2 tau)"
    grouplabels[6]="4-f (>2 tau)"

    plotfilename=figdir
    if saveFitHistos and not plotHistos:
        plotfilename= plotfilename+'FIT'
    plotfilename= plotfilename+'tautau_'+model+pol+'.pdf'

    #print orderedClasses
    #for cl in orderedClasses:
    #    print cl, rootfiles[cl].keys()

    legend = 0

    cc2 = ROOT.TCanvas("cc2","cc2",500,500)
    cc = ROOT.TCanvas()
    cc.Print(plotfilename+'[')

    tt = ROOT.TText()
    tt.SetTextColor(2)

    nplots=0

    plotstokeep=[]
    xaxistitle={}

    if plotHistos:
        plotstokeep.append("sample0_visjetAcolinearity")
        xaxistitle["sample0_visjetAcolinearity"]='acolinearity [rad]'
        plotstokeep.append("sample0_visjetAcoplanarity")
        xaxistitle["sample0_visjetAcoplanarity"]='acoplanarity [rad]'
        plotstokeep.append("sample0_outsideEnergy")
        xaxistitle["sample0_outsideEnergy"]='energy outside cones [GeV]'
        plotstokeep.append("sample0_outsidePt")
        xaxistitle["sample0_outsidePt"]='p_{T} outside cones [GeV/c]'
        plotstokeep.append("sample0_minSeedEn")
        xaxistitle["sample0_minSeedEn"]='smaller seed energy [GeV]'
        plotstokeep.append("sample0_maxSeedEn")
        xaxistitle["sample0_maxSeedEn"]='larger seed energy [GeV]'

        plotstokeep.append("sample0_ooconeMaxGammaEn")
        xaxistitle["sample0_ooconeMaxGammaEn"]='max. energy photon PFO ouside cone [GeV]'

        plotstokeep.append("sample0_ooconeMaxPFOEn")
        xaxistitle["sample0_ooconeMaxPFOEn"]='max. energy PFO ouside cone [GeV]'

        plotstokeep.append("sample0__ALL_SELcone_ngammapfo")
        xaxistitle["sample0__ALL_SELcone_ngammapfo"]='number of photons in cone [GeV]'

        plotstokeep.append("sample0__ALL_SELconeTRIM_visMass")
        xaxistitle["sample0__ALL_SELconeTRIM_visMass"]='visible mass of cone [GeV/c~{2}]'

        plotstokeep.append("sample0__ALL_SELconeTRIM_neutralvisMass")
        xaxistitle["sample0__ALL_SELconeTRIM_neutralvisMass"]='visible neutral mass of cone [GeV/c~{2}]'

        plotstokeep.append("sample0_SEL_coneMass")
        xaxistitle["sample0_SEL_coneMass"]='visible mass of cone [GeV/c~{2}]'

        plotstokeep.append("sample0_SEL_rec_pi_pol")
        xaxistitle["sample0_SEL_rec_pi_pol"]='polarimeter reco #tau -> #pi decays'

        plotstokeep.append("sample0_SEL_rec_pi_coneMass")
        xaxistitle["sample0_SEL_rec_pi_coneMass"]='visible mass reco #tau -> #pi decays'

        plotstokeep.append("sample0_SEL_rec_rho_pol")
        xaxistitle["sample0_SEL_rec_rho_pol"]='polarimeter reco #tau -> #rho decays'

        plotstokeep.append("sample0_SEL_rec_rho_coneMass")
        xaxistitle["sample0_SEL_rec_rho_coneMass"]='visible mass reco #tau -> #rho decays'

        plotstokeep.append("sample0_SEL_rec_a1p_coneMass")
        xaxistitle["sample0_SEL_rec_a1p_pol"]='visible mass reco #tau -> a1p decays'

        plotstokeep.append("sample0_logseedClusterWidth1")
        xaxistitle["sample0_logseedClusterWidth1"]="log_{10} ( cluster ellipse eval_{min} [mm] )"

        plotstokeep.append("sample0_logseedClusterLength")
        xaxistitle["sample0_logseedClusterLength"]="log_{10} ( cluster ellipse eval_{max} [mm] )"

    patternsNotToPlot=[]
    patternsNotToPlot.append('sample0__OTHER_')
    patternsNotToPlot.append('sample0__M_')
    patternsNotToPlot.append('sample0__E_')


    for histoname in histonames2:
        ignore=False
        for patt in patternsNotToPlot:
            if patt in histoname:
                ignore=True
                break
        if ignore:
            continue
        groupedHistos=[]

        #print orderedClasses

        if 'sample0' in histoname:

            # print histoname

            for cl in orderedClasses:
                # print cl, rootfiles[cl].keys()
                for rk in rootfiles[cl].keys():
                    rootfile = rootfiles[cl][rk]
                    hh = rootfile.Get(histoname)

                    if len( groupedHistos ) == 0 :
                        for ig in range(7):
                            newname=histoname.replace('sample0','group'+str(ig) ) 
                            hf = hh.Clone( newname )
                            hf.SetName( newname )
                            # hf.Clear()
                            hf.Reset()
                            hf.SetLineColor(cols[ig])
                            hf.SetFillColor(cols[ig])
                            hf.SetFillStyle(1001)
                            groupedHistos.append( hf )

                    if legend==0:
                        legend = ROOT.TLegend(0.75, 0.65, 0.95, 0.95)
                        for igh in range(7):
                            legend.AddEntry( groupedHistos[igh], grouplabels[igh], "f" )
            
                    group=-1
                    for gr in grouping.keys():
                        if gr in cl:
                            group=gr
                            break

                    for icl in range(nclass):
                        hh = rootfile.Get(histoname.replace('sample0','sample'+str(icl)) )
                        hh.Scale(procweight[rk])
                        newgroup = grouping[group][icl]
                        groupedHistos[newgroup].Add( hh )

            hstack = ROOT.THStack("STACK"+histoname, histoname)
            hstack2 = ROOT.THStack("STACK2"+histoname, histoname)
            for gh in range(7):
                hstack.Add(groupedHistos[gh])
            for gh in range(6,-1,-1):
                hstack2.Add(groupedHistos[gh])


            if saveFitHistos and 'SEL_rec_' in histoname:
                print 'trying to write histos to file...', histoname
                fithistos.cd()
                for gh in groupedHistos:
                    gh.Write()
                hstack.Write()
                hstack2.Write()


            cc.Clear()
            cc.Divide(4,3)
            for gh in range(7):
                cc.cd(gh+1)
                groupedHistos[gh].Draw("hist")
            cc.cd(8).SetLogy()
            hstack.SetMinimum(1)
            hstack.Draw("hist")
            legend.Draw()
            cc.cd(9)
            hstack.Draw("hist")
            legend.Draw()

            cc.cd(10).SetLogy()
            hstack2.SetMinimum(1)
            hstack2.Draw("hist")
            legend.Draw()
            cc.cd(11)
            hstack2.Draw("hist")
            legend.Draw()

            cc.Print(plotfilename, 'Title:'+histoname)

            if histoname in plotstokeep:

                xtit=histoname
                if histoname in xaxistitle.keys():
                    xtit=xaxistitle[histoname]

                hstack.GetXaxis().SetTitle(xtit)
                hstack2.GetXaxis().SetTitle(xtit)

                cc2.Clear()
                cc2.SetLogy(0)
                hstack.SetTitle('')
                hstack.SetMinimum(0)
                hstack.Draw("hist")
                legend.Draw()
                ttf.DrawTextNDC(labPosX, labPosY, ILDlab)
                ttf2.DrawTextNDC(lab2PosX, lab2PosY, beamLab)
                cc2.Print(figdir+'/'+histoname+'_a.eps')
                cc2.Print(figdir+'/'+histoname+'_a.C')

                cc2.Clear()
                cc2.SetLogy(1)
                hstack.SetMinimum(1)
                hstack.Draw("hist")
                legend.Draw()
                ttf.DrawTextNDC(labPosX, labPosY, ILDlab)
                ttf2.DrawTextNDC(lab2PosX, lab2PosY, beamLab)
                cc2.Print(figdir+'/'+histoname+'_b.eps')
                cc2.Print(figdir+'/'+histoname+'_b.C')

                cc2.Clear()
                cc2.SetLogy(0)
                hstack2.SetTitle('')
                hstack2.SetMinimum(0)
                hstack2.Draw("hist")
                legend.Draw()
                ttf.DrawTextNDC(labPosX, labPosY, ILDlab)
                ttf2.DrawTextNDC(lab2PosX, lab2PosY, beamLab)
                cc2.Print(figdir+'/'+histoname+'_c.eps')
                cc2.Print(figdir+'/'+histoname+'_c.C')

                cc2.Clear()
                cc2.SetLogy(1)
                hstack2.SetMinimum(1)
                hstack2.Draw("hist")
                legend.Draw()
                ttf.DrawTextNDC(labPosX, labPosY, ILDlab)
                ttf2.DrawTextNDC(lab2PosX, lab2PosY, beamLab)
                cc2.Print(figdir+'/'+histoname+'_d.eps')
                cc2.Print(figdir+'/'+histoname+'_d.C')

    cc.Print(plotfilename+']')
    print 'closed pdf file', plotfilename

    for cl in rootfiles.keys():
        tothist=0
        for rk in rootfiles[cl].keys():
            rootfiles[cl][rk].Close()

if saveFitHistos:
    fithistos.Close()
                        
