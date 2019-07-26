import os
import subprocess

topindir='/group/ilc/soft/samples/mc-opt-3/ild/dst-merged/500-TDR_ws/'


models=['ILD_l5_o1_v02', 'ILD_s5_o1_v02']

samples=['2f_Z_leptonic']

# samples=['2f_Z_bhabhag','2f_Z_hadronic','2f_Z_leptonic','4f_WW_hadronic','4f_WW_leptonic','4f_WW_semileptonic','4f_ZZWWMix_hadronic','4f_ZZWWMix_leptonic','4f_ZZ_hadronic','4f_ZZ_leptonic','4f_ZZ_semileptonic','4f_singleW_leptonic','4f_singleW_semileptonic','4f_singleZee_leptonic','4f_singleZee_semileptonic','4f_singleZnunu_leptonic','4f_singleZnunu_semileptonic','4f_singleZsingleWMix_leptonic']
# ,'5f','6f_eeWW','6f_llWW','6f_ttbar','6f_vvWW','6f_xxWW','6f_xxxxZ','6f_yyyyZ','aa_4f','higgs_ffh']

maxevt=0
# maxevt=10000000


# samples=['2f_Z_bhabhag','2f_Z_hadronic','2f_Z_leptonic']
# pols=['eL.pR']

for model in models:
    modeldir=model+'/v02-00-01/'

    for samp in samples:
        print model, samp

        indir=topindir+samp+'/'+modeldir

        procsPols={}
        inff={}
        infiles=os.listdir(indir)
        for infile in infiles:
            processid = infile.split('.')[4]
            pol=infile.split('.')[6]+'.'+infile.split('.')[7]
            if processid not in procsPols.keys():
                procsPols[processid]={}
            if pol not in procsPols[processid]:
                procsPols[processid][pol]=[]
            procsPols[processid][pol].append(infile)

        for kk in  procsPols.keys():
            print kk, procsPols[kk].keys()

        for proc in procsPols.keys():

            for pol in procsPols[proc]:
        
                jobname=model+'_'+samp+'_'+proc+'_'+pol+'_oneHadSS'

                if maxevt==0:
                    jobname='ALLEVT_'+jobname

                steername='steer_'+jobname+'.xml'
                logname='log_'+jobname+'.txt'

                steerfile = open(steername, 'w')
                steerfile.write('<marlin> \n') 
                steerfile.write('  <execute> \n') 
                steerfile.write('    <processor name="InitDD4hep"/> \n') 
                steerfile.write('    <processor name="MyIsolatedLeptonTaggingProcessor"/> \n')
                steerfile.write('    <processor name="myTauTauFinder"/> \n') 
            # steerfile.write('    <processor name="myAnaTauTau"/> \n') 
                steerfile.write('  </execute> \n') 
                steerfile.write('  <global> \n')
                steerfile.write('    <parameter name="LCIOInputFiles"> \n')

                inff=procsPols[proc][pol]
                inff.sort()
                for infff in inff:
                    steerfile.write(indir+'/'+infff+'\n')

                steerfile.write('    </parameter>\n')
                steerfile.write('    <parameter name="MaxRecordNumber" value="'+str(maxevt)+'"/>\n')
                steerfile.write('    <parameter name="SkipNEvents" value="0"/>\n')
                steerfile.write('    <parameter name="SupressCheck" value="false"/>\n')
                steerfile.write('    <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT"> MESSAGE </parameter>\n')
                steerfile.write('  </global>\n')
                steerfile.write('  <processor name="InitDD4hep" type="InitializeDD4hep">\n')
                steerfile.write('    <parameter name="DD4hepXMLFile" type="string"> $lcgeo_DIR/ILD/compact/'+model+'/'+model+'.xml </parameter>\n')
                steerfile.write('  </processor>\n')
                steerfile.write('  <processor name="myTauTauFinder" type="danielKeitaTauFinderProcessor">\n')
                steerfile.write('     <parameter name="outputFilename"> tauFind_'+jobname+' </parameter>\n')
                steerfile.write('  </processor>\n')
                steerfile.write('  <processor name="myAnaTauTau"  type="analyseTauProcessor">\n')
                steerfile.write('     <parameter name="outputFilename"> tauAna_'+jobname+'.root </parameter>\n')
                steerfile.write('  </processor>\n')

                steerfile.write('  <processor name="MyIsolatedLeptonTaggingProcessor" type="IsolatedLeptonTaggingProcessor">\n')
                steerfile.write('   <parameter name="CosConeLarge" type="float">0.95 </parameter>\n')
                steerfile.write('   <parameter name="CosConeSmall" type="float">0.98 </parameter>\n')
                steerfile.write('   <parameter name="CutOnTheISOElectronMVA" type="float">0.8 </parameter>\n')
                steerfile.write('   <parameter name="CutOnTheISOMuonMVA" type="float">0.8 </parameter>\n')
                steerfile.write('   <parameter name="DirOfISOElectronWeights" type="string">isolated_electron_weights </parameter>\n')
                steerfile.write('   <parameter name="DirOfISOMuonWeights" type="string">isolated_muon_weights_woYoke </parameter>\n')
                steerfile.write('   <parameter name="InputPandoraPFOsCollection" type="string" lcioInType="ReconstructedParticle">PandoraPFOs </parameter>\n')
                steerfile.write('   <parameter name="InputPrimaryVertexCollection" type="string" lcioInType="ReconstructedParticle">PrimaryVertex </parameter>\n')
                steerfile.write('   <parameter name="IsSelectingOneIsoLep" type="bool"> false </parameter>\n')
                steerfile.write('   <parameter name="MaxD0SigForElectron" type="float">50 </parameter>\n')
                steerfile.write('   <parameter name="MaxD0SigForMuon" type="float">20 </parameter>\n')
                steerfile.write('   <parameter name="MaxEOverPForElectron" type="float">1.3 </parameter>\n')
                steerfile.write('   <parameter name="MaxEOverPForMuon" type="float">0.3 </parameter>\n')
                steerfile.write('   <parameter name="MaxZ0SigForElectron" type="float">50 </parameter>\n')
                steerfile.write('   <parameter name="MaxZ0SigForMuon" type="float">20 </parameter>\n')
                steerfile.write('   <parameter name="MinEOverPForElectron" type="float">0.5 </parameter>\n')
                steerfile.write('   <parameter name="MinEecalOverTotEForElectron" type="float">0.9 </parameter>\n')
                steerfile.write('   <parameter name="MinPForElectron" type="float">5 </parameter>\n')
                steerfile.write('   <parameter name="MinPForMuon" type="float">5 </parameter>\n')
                steerfile.write('   <parameter name="OutputIsoLeptonsCollection" type="string" lcioOutType="ReconstructedParticle">ISOLeptons </parameter>\n')
                steerfile.write('   <parameter name="OutputPFOsWithoutIsoLepCollection" type="string" lcioOutType="ReconstructedParticle">PandoraPFOsWithoutIsoLep </parameter>\n')
                steerfile.write('  </processor>\n')

                steerfile.write('</marlin>\n')

                steerfile.close()

                command = 'bsub -q s "Marlin ' + steername + ' &> ' + logname + '"'

                print command

                ret=subprocess.call(command, shell=True)

