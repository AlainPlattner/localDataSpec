function spec_expected = expectedSpecML(specML,Lfull,Ltap)

lsFull = 0:Lfull;
% First transfom from ML to powerspectral density
specPD = specML(:)./(lsFull(:)+1)./(2*lsFull(:)+1);
% LocalizeIt
specPD_loc = localizeSpec(specPD,Ltap)';
% Back to ML
spec_expected = specPD_loc(:).*(lsFull(:)+1).*(2*lsFull(:)+1);
