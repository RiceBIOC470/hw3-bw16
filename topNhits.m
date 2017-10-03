function [ accession_numbers ] = topNhits( accession,N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
gb_data=getgenbank(accession);
seq=gb_data.Sequence;
[requestID,requestTime]=blastncbi(seq,'blastn');
blast_data=getblast(requestID,'WaitTime',requestTime);
name={blast_data.Hits.Name};
name_N=name(1:N);

for i=1:N
    name_split=strsplit(name_N{i},'|');
    accession_numbers(i)=name_split(4);
end
end

