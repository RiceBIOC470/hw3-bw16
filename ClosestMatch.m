function [ h,nh ] = ClosestMatch(accession)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
accessionNumbers=topNhits(accession,50);

subject=getgenbank(accession);
humanList=cell(50,1);
nonhumanList=cell(50,1);
for ii=1:length(accessionNumbers)
info=getgenbank(accessionNumbers{ii});
species=info.Source;
a=strfind(species,'Homo sapiens (human)');
if a==1
    humanList{ii}=accessionNumbers{ii};
else
    nonhumanList{ii}=accessionNumbers{ii};
    
end
end
highestHscore=0;
if isempty(humanList)==1
    disp('No human homologue found.');
else 
    humanList = humanList(~cellfun('isempty',humanList))
 for jj=1:length(humanList)
     queryh=getgenbank(humanList{jj});
     [score] = swalign(subject.Sequence,queryh.Sequence,'Alphabet','nt','Showscore',false);
     if score>highestHscore
         highestHscore=score;
         h=humanList{jj};
     end
 end
end

highestNHscore=0
if isempty(nonhumanList)==1
    disp('No non-human homologue found');
else
    nonhumanList = nonhumanList(~cellfun('isempty',nonhumanList))
    for xx=1:length(nonhumanList)
     querynh=getgenbank(nonhumanList{xx});
     [score] = swalign(subject.Sequence,querynh.Sequence,'Alphabet','nt','Showscore',false);
     if score>highestNHscore
         highestNHscore=score;
         nh=nonhumanList{xx};
     end
 end
end
hgenbank=getgenbank(h);
nhgenbank=getgenbank(nh);
human=[{hgenbank.Accession},{hgenbank.Definition}]
nonhuman=[{nhgenbank.Accession},{nhgenbank.Definition}]
end
    
    

