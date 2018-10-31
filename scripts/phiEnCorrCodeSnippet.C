float getBarrelFoldedPhi( float phi) {
  if ( phi<0 ) phi+=2*acos(-1);
  int sector = phi/(acos(-1)/4);
  float foldedPhi = phi - sector*acos(-1)/4.;
  return foldedPhi;
}

float phiCorrectedEnergy( float inEn, float inPhi ) {  

  // function to correct phi-dependence of photon PFO energies IN BARREL REGION ONLY
  // inEn  = PFO energy [GeV]
  // inPhi = PFO azimuthal angle [rad]

  float logen = log10(inEn);

  float foldedPhi = getBarrelFoldedPhi( inPhi );

  float par1 = 0.983607   + logen*0.018598;  // overall scale
  float par2 = 0.414583   + logen*0.0126418; // position in folded phi
  float par3 = 0.0995     - logen*0.0033125; // depth of dip
  float par4 = 0.00819825 + logen*0.00238604; // width of dip (neg side)
  float par5 = 0.0344557  + logen*0.00350749; // width of dip (pos side)

  float dipdist = foldedPhi - par2;
  float sig = dipdist < 0 ? par4 : par5;

  float corfn = par1 - par3*exp( -pow(dipdist,2)/(2*pow(sig,2) ) );

  return inEn/corfn;
}
