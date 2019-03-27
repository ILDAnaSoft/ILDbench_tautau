import os
import subprocess
import ROOT


finLR = ROOT.TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250106_eL.pR.root", "read");
finRL = ROOT.TFile("tauFind_ILD_l5_o1_v02_2f_Z_leptonic_I250108_eR.pL.root", "read");


#pol='unpol'
pol='eLpR'
# pol='eRpL'

totlumi=4000 # in fb-1

#--------------------------

plotHistos=True

#--------------------------
elePol=0.
posPol=0.
lumi=0
if pol=='eLpR':
    elePol=-0.8
    posPol=0.3
    lumi=0.4*totlumi
elif pol=='eRpL':
    elePol=0.8
    posPol=-0.3
    lumi=0.4*totlumi
elif pol=='eRpR':
    elePol=0.8
    posPol=0.3
    lumi=0.1*totlumi
elif pol=='eLpL':
    elePol=-0.8
    posPol=-0.3
    lumi=0.1*totlumi
else:
    pol='unpol'

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

ntoy=500



print '**********************************************'
print '    elePol=', elePol, 'posPol=', posPol, 'lumi=', lumi
print '**********************************************'

cc=ROOT.TCanvas()
trand=ROOT.TRandom()


cc.Print('mcana.pdf[')


hLRorig=finLR.Get("mctautau_ecom_tauMinusCosth")
hRLorig=finRL.Get("mctautau_ecom_tauMinusCosth")

cc.Divide(5,3)
cc.cd(1)
hLRorig.Draw("zcol")

cc.cd(2)
hRLorig.Draw("zcol")

# fitfn=ROOT.TF1("costhfn","[0]*(1 + [1]*x + [2]*x*x)/(2*(1.+[2]/3.))", -1, 1 )
fitfn=ROOT.TF1("costhfn","[0]*( [1]*(1+x)*(1+x)  + (1-[1])*(1-x)*(1-x) )", -1, 1 )

pp=[]

for ipol in range(2):
    
    hLR=hLRorig.Clone('LRpol'+str(ipol))
    hRL=hRLorig.Clone('RLpol'+str(ipol))
    polstr=''

    if ipol==0:
        polstr='BeamPol_eL80pR30'
        hLR.Scale( wt_eL*wt_pR )
        hRL.Scale( wt_eR*wt_pL )
    else:
        polstr='BeamPol_eR80pL30'
        hLR.Scale( wt_eR*wt_pL )
        hRL.Scale( wt_eL*wt_pR )

    cc.cd(5*(ipol+1)+1)

    ibin1=hLR.GetXaxis().FindBin(480.)
    projLR = hLR.ProjectionY( "_LR_"+polstr, hLR.GetXaxis().FindBin(480), hLR.GetNbinsX() )
    projLR.Fit("costhfn","L")
    projLR.SetMinimum(0)
    projLR.Draw()

    print fitfn.Integral(-1,1)/(projLR.GetBinWidth(2)),  # divide by bin width
    print projLR.Integral()

    cc.cd(5*(ipol+1)+2)
    projRL = hRL.ProjectionY( "_RL_"+polstr, hRL.GetXaxis().FindBin(480), hRL.GetNbinsX() )
    projRL.Fit("costhfn","L")
    projRL.SetMinimum(0)
    projRL.Draw()

    cc.cd(5*(ipol+1)+3)
    projSum = projLR.Clone("Total_"+polstr)
    projSum.Add(projRL)
    projSum.Fit("costhfn","L")
    projSum.SetMinimum(0)
    projSum.Draw()
    pp.append(projSum)


    bigFitPar0 = projSum.GetFunction("costhfn").GetParameter(0)
    bigFitPar1 = projSum.GetFunction("costhfn").GetParameter(1)

    hFitErr = ROOT.TH1F("fiterr","fiterr",100,0,0.01)
    hFitPull = ROOT.TH1F("fiterrpull","fiterrpull",100,-5,5)

    htoy = projSum.Clone("toy")
    for i in range(ntoy):
        htoy.Reset()
        fitfn.SetParameter(0, bigFitPar0)
        fitfn.SetParameter(1, bigFitPar1)
        meanevt = projSum.GetFunction("costhfn").Integral(-1,1)/projSum.GetBinWidth(3)
        nevt = trand.Poisson( meanevt )
        htoy.FillRandom( "costhfn", nevt )
        htoy.Fit( "costhfn" ,"qL"  )
        fpar=htoy.GetFunction("costhfn").GetParameter(1)
        fparErr=htoy.GetFunction("costhfn").GetParError(1)
        hFitErr.Fill(fparErr)
        hFitPull.Fill( (fpar-bigFitPar1)/fparErr )

    cc.cd(5*(ipol+1)+4)
    hFitErr.Draw()

    cc.cd(5*(ipol+1)+5)
    hFitPull.Draw()

    pp.append(htoy)
    pp.append(hFitErr)
    pp.append(hFitPull)






cc.Print('mcana.pdf')

cc.Clear()

hLRorig2d=finLR.Get("sample0_pirho_mcPolarExactPlusMinus")
hRLorig2d=finRL.Get("sample0_pirho_mcPolarExactPlusMinus")



hLRorig=finLR.Get("sample0_pi_mcPolarExact")
hRLorig=finRL.Get("sample0_pi_mcPolarExact")

cc.Divide(5,5)
cc.cd(1)
hLRorig2d.Draw("zcol")

cc.cd(2)
hRLorig2d.Draw("zcol")

cc.cd(3)
hLRorig.SetMinimum(0)
hLRorig.Draw()

cc.cd(4)
hRLorig.SetMinimum(0)
hRLorig.Draw()

ic=6

fitpol = ROOT.TF1("fpol","[0]*( [1]*(1+x) + (1.-[1])*(1-x) )", -1, 1)


for ipr in range(2):

    if ipr==0:
        hLRorig=finLR.Get("sample0_pi_mcPolarExact")
        hRLorig=finRL.Get("sample0_pi_mcPolarExact")
    else:
        hLRorig=finLR.Get("sample0_rho_mcPolarExact")
        hRLorig=finRL.Get("sample0_rho_mcPolarExact")

    for ipol in range(2):
    
        hLR=hLRorig.Clone('LRpol'+str(ipol))
        hRL=hRLorig.Clone('RLpol'+str(ipol))
        polstr=''

        if ipol==0:
            polstr='BeamPol_eL80pR30'
            hLR.Scale( wt_eL*wt_pR )
            hRL.Scale( wt_eR*wt_pL )
        else:
            polstr='BeamPol_eR80pL30'
            hLR.Scale( wt_eR*wt_pL )
            hRL.Scale( wt_eL*wt_pR )

        cc.cd( ic )

        hLR.SetMinimum(0)
        hLR.Draw()

        ic=ic+1
        cc.cd( ic )

        hRL.SetMinimum(0)
        hRL.Draw()

        ic=ic+1
        cc.cd( ic )

        hSum = hLR.Clone("Total_"+polstr)
        hSum.Add(hRL)
        hSum.SetMinimum(0)
        hSum.Fit("fpol","L")
        hSum.Draw()

        ic=ic+1
        cc.cd( ic )


        pp.append(hLR)
        pp.append(hRL)
        pp.append(hSum)

        bigFitPar0 = hSum.GetFunction("fpol").GetParameter(0)
        bigFitPar1 = hSum.GetFunction("fpol").GetParameter(1)

        hFitErr = ROOT.TH1F("fiterr","fiterr",100,0,0.03)
        hFitPull = ROOT.TH1F("fiterrpull","fiterrpull",100,-5,5)

        htoy = hSum.Clone("toy")
        for i in range(ntoy):
            htoy.Reset()
            fitpol.SetParameter(0, bigFitPar0)
            fitpol.SetParameter(1, bigFitPar1)
            meanevt = hSum.GetFunction("fpol").Integral(-1,1)/hSum.GetBinWidth(3)
            nevt = trand.Poisson( meanevt )
            htoy.FillRandom( "fpol", nevt )
            htoy.Fit( "fpol" ,"qL"  )
            fpar=htoy.GetFunction("fpol").GetParameter(1)
            fparErr=htoy.GetFunction("fpol").GetParError(1)
            hFitErr.Fill(fparErr)
            hFitPull.Fill( (fpar-bigFitPar1)/fparErr )

        hFitErr.Draw()

        ic=ic+1
        cc.cd( ic )
        hFitPull.Draw()

        ic=ic+1
        cc.cd( ic )

        pp.append(htoy)
        pp.append(hFitErr)
        pp.append(hFitPull)




cc.Print('mcana.pdf')



cc.Print('mcana.pdf]')
