function V = parseVCF(filename)
verbose = true;
%% function V = parseVCF(filename)
% #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
% 20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.2
V = [];
tic;
try 
    fdin = fopen(filename);
catch fex
    disp(['unable to open ' filename]);
    return;
end

nhdr = 0;
while 1
    hdrstr = fgets(fdin);
%     if strcmp('chr',hdrstr(1:3)); % warning: chr not std on all refs
%         frewind(fdbed);
%         break;
%     end;
    if ~isempty(hdrstr) && hdrstr(1)~='#'  % start data - not comment or blank
        examplerec= hdrstr;
        frewind(fdin);
        break;
    end
    nhdr = nhdr+1;
end
V.hdr = cell(nhdr,1);
for i = 1:nhdr
    V.hdr{i} = fgets(fdin);
    if verbose
        fprintf(V.hdr{i});
    end
end
% determine number of samples in file
% #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
f = strsplit(V.hdr{end},'\t');
f{end} = f{end}(1:end-1);
nsamp = numel(f)-9;
fmt = ['%s %d32 %s %s %s %f32 %s %s %s' repmat(' %s',nsamp) ' %*[^\n]'];
D = textscan(fdin, fmt, 'delimiter','\t');
toc;
fclose(fdin);
V.R(:,3) = D{2}; D{2} = [];
V.R(:,2) = V.R(:,3);
[V.segNames, ~, V.R(:,1)] = unique(D{1}); D{1} = [];

%% compression schemes (may not all be worthwhile)
V.ID = D{3}; D{3} = [];
[V.uREF, ~, V.jREF] = unique(D{4}); D{4} = []; V.jREF = int32(V.jREF); 
[V.uALT, ~, V.jALT] = unique(D{5}); D{5} = []; V.jALT = int32(V.jALT);
V.QUAL = D{6};
[V.uFILTER, ~, V.jFILTER] = unique(D{7}); D{7} = []; V.jFILTER = int32(V.jFILTER);
% Info is a potential memory hog and should be optional
% [V.uINFO, ~, V.jINFO] = unique(D{8}); D{8} = []; V.jINFO = int32(V.jINFO);
V.INFO = D{8};
[V.uFORMAT, ~, V.jFORMAT] = unique(D{9}); D{9} = []; V.jFORMAT = int32(V.jFORMAT);
for i = 10:numel(f)
    V.(['samp' f{i}]) = D{i};
end
clear D;
toc;






