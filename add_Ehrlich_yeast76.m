function [yeast76ehrlich] = add_Ehrlich_yeast76(yeast76)

yeast76ehrlich = yeast76;
%add Ehrlich pathway catabolism for L-tyrosine
%YGL202W ARO8
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'L-tyrosine transaminase (L-glu)',{yeast76.mets{850},yeast76.mets{124},'4-hydroxyphenylpyruvate [cytosolic]',yeast76.mets{792}},[-1 -1 1 1],1,-1000,1000,0,[],{'YGL202W'},{'YGL202W'},{'YGL202W'},'true');
%YHR137W ARO9
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'L-tyrosine transaminase (L-ala)',{yeast76.mets{850},yeast76.mets{1088},'4-hydroxyphenylpyruvate [cytosolic]',yeast76.mets{756}},[-1 -1 1 1],1,-1000,1000,0,[],{'YHR137W'},{'YHR137W'},{'YHR137W'},'true');

%(YGR087C OR YLR044C OR YLR134W) PDC6/1/5
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'4-hydroxyphenylpyruvate decarboxylase',{'4-hydroxyphenylpyruvate [cytosolic]',yeast76.mets{605},yeast76.mets{343},'4-hydroxyphenylacetaldehyde [cytosolic]'},[-1 -1 1 1],0,0,1000,0,[],{'YGR087C OR YLR044C OR YLR134W'},{'YGR087C','YLR044C','YLR134W'},{'YGR087C','YLR044C','YLR134W'},'true');

%YPL061W ALD6
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'aldehyde dehydrogenase (4-hydroxyphenylacetate)',{'4-hydroxyphenylacetaldehyde [cytosolic]',yeast76.mets{614},yeast76.mets{961},'4-hydroxyphenylacetate [cytosolic]',yeast76.mets{605},yeast76.mets{965}},[-1 -1 -1 1 1 1],1,-1000,1000,0,[],{'YPL061W'},{'YPL061W'},{'YPL061W'},'true');

%ADH5, SFA1, ADH1 (YBR145W OR YDL168W OR YOL086C) 
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'aldehyde dehydrogenase (tyrosol)',{'4-hydroxyphenylacetaldehyde [cytosolic]',yeast76.mets{952},'tyrosol [cytosolic]',yeast76.mets{605},yeast76.mets{957}},[-1 1 1 -1 -1],1,-1000,1000,0,[],{'YBR145W OR YDL168W OR YOL086C'},{'YBR145W','YDL168W','YOL086C'},{'YBR145W','YDL168W','YOL086C'},'true');

[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'tyrosol transport',{'tyrosol [extracellular]','tyrosol [cytosolic]'},[1 -1],1,-1000,1000);
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'tyrosol exchange',{'tyrosol [extracellular]'},[-1],0,0,1000);

[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'4-hydroxyphenylacetate transport',{'4-hydroxyphenylacetate [extracellular]','4-hydroxyphenylacetate [cytosolic]'},[1 -1],1,-1000,1000);
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'4-hydroxyphenylacetate exchange',{'4-hydroxyphenylacetate [extracellular]'},[-1],0,0,1000);

[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'4-hydroxyphenylpyruvate transport',{'4-hydroxyphenylpyruvate [extracellular]','4-hydroxyphenylpyruvate [cytosolic]'},[1 -1],1,-1000,1000);
[yeast76ehrlich,rxnIDexists] = addReaction(yeast76ehrlich,'4-hydroxyphenylpyruvate exchange',{'4-hydroxyphenylpyruvate [extracellular]'},[-1],0,0,1000);

end

