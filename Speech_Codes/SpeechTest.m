% American Sign Language Detection-Speech
% Performs testing on the testing samples of the speech database.
clc
clear all
close all

%% Add training and test directories to path
addpath('./Training/')
cAlpha = [{'A'},{'B'},{'C'},{'Five'},{'Point'},{'V'}];%No of alphabet used

%% Load Training set
load Fan_all.mat
load Neeraj_all.mat
load Zhang_all.mat

%% Set variables
NoOfSamples = 20;
NoOfTestSamples = 5;            
Users = {'Fan','Neeraj','Zhang'};
Letters = {'A','B','C','Five','Point','V'};

%MFCC parameters
OverlapSize = 0.5;
MFCCNo = 45;
NoOfWindows = 25;
NoOfFilters = floor(MFCCNo/NoOfWindows+1);
mfcc = zeros(size(Letters,2),size(Users,2)*NoOfSamples,MFCCNo);

%% MFCC Training for all letters
for ii = 1:size(Letters,2);    
    for jj = 1:size(Users,2);
        for kk = 1:NoOfSamples            
            file_name = strcat(Users(jj),'_',Letters(ii),int2str(kk));
            Samples = eval(char(file_name));
            zz = find(Samples) < max(Samples/3);%Threshold speech regions
            Samples(zz) = 0;
            zz = find(Samples);
            Speech_Region = Samples(zz)/norm(Samples(zz));             
            WindowSize = floor((size(Speech_Region,1))/(NoOfWindows+1));
            ww = 0;
            for ll = 0:OverlapSize:(NoOfWindows-1)/2
                bb = Speech_Region(floor(ll*WindowSize)+1:floor(ll*WindowSize)+WindowSize).*hamming(WindowSize);
                fb = fft(bb);                
                mb = 2595 * log10(1 + fb./700);                
                mfout = dct(log(abs(mb)),NoOfFilters);                                    
                mfcc(ii,kk,ww*NoOfFilters+1:ww*NoOfFilters+NoOfFilters) = mfout;                                                        
                ww = ww + 1;
            end            
        end
    end
end


%% Perform Gaussian Modelling for MFCC
Windows = size(mfcc,3);
tempStorage = zeros(size(Users,2)*NoOfSamples,Windows);
tempStorage(:,:) = mfcc(1,:,:);
obj_A = gmdistribution.fit(tempStorage,1,'Regularize',0.01);
tempStorage(:,:) = mfcc(2,:,:);
obj_B = gmdistribution.fit(tempStorage,1,'Regularize',0.01);
tempStorage(:,:) = mfcc(3,:,:);
obj_C = gmdistribution.fit(tempStorage,1,'Regularize',0.01);
tempStorage(:,:) = mfcc(4,:,:);
obj_Five = gmdistribution.fit(tempStorage,1,'Regularize',0.01);
tempStorage(:,:) = mfcc(5,:,:);
obj_Point = gmdistribution.fit(tempStorage,1,'Regularize',0.01);
tempStorage(:,:) = mfcc(6,:,:);
obj_V = gmdistribution.fit(tempStorage,1,'Regularize',0.01);

% LPC parameters
NoOfLPCFilters = 50;
lpccoeff = zeros(size(Letters,2),size(Users,2)*NoOfSamples,NoOfLPCFilters+1);

%% LPC Training for all letters
% Training for all letters
for ii = 1:size(Letters,2);
    ll = 1;
    for jj = 1:size(Users,2);
        for kk = 1:NoOfSamples            
            file_name = strcat(Users(jj),'_',Letters(ii),int2str(kk));
            Samples = eval(char(file_name));
            zz = find(Samples) < max(Samples/3);%Threshold speech regions
            Samples(zz) = 0;
            zz = find(Samples);
            Speech_Region = Samples(zz)/norm(Samples(zz));            
            lpccoeff(ii,ll,:) = lpc(Speech_Region,NoOfLPCFilters);            
            ll = ll + 1;
        end
    end
end

%% Perform Gaussian Modelling for LPC
tempStorage = zeros(size(Users,2)*NoOfSamples,NoOfLPCFilters);
tempStorage(:,:) = lpccoeff(1,:,2:end);
obj_A2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(2,:,2:end);
obj_B2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(3,:,2:end);
obj_C2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(4,:,2:end);
obj_Five2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(5,:,2:end);
obj_Point2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(6,:,2:end);
obj_V2 = gmdistribution.fit(tempStorage,1);

%% Perform Gaussian Modelling for LPC
tempStorage = zeros(size(Users,2)*NoOfSamples,NoOfLPCFilters);
tempStorage(:,:) = lpccoeff(1,:,2:end);
obj_A2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(2,:,2:end);
obj_B2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(3,:,2:end);
obj_C2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(4,:,2:end);
obj_Five2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(5,:,2:end);
obj_Point2 = gmdistribution.fit(tempStorage,1);
tempStorage(:,:) = lpccoeff(6,:,2:end);
obj_V2 = gmdistribution.fit(tempStorage,1);

LPCCorrect_Test = zeros(size(Users,2),size(Letters,2));
MFCCCorrect_Test = zeros(size(Users,2),size(Letters,2));
for ii = 1:size(Letters,2);
    for jj = 1:size(Users,2);
        for kk = 1:NoOfTestSamples            
			%% Extract MFCC for test data
			InputIn = strcat(Users(jj),'_',Letters(ii),int2str(kk));
			Samples = eval(char(InputIn));
			zz = find(Samples) < max(Samples/3);%Threshold speech regions
			Samples(zz) = 0;
			zz = find(Samples);
			Speech_Region = Samples(zz);  
			mfcc_test = zeros(1,MFCCNo);
			WindowSize = floor((size(Speech_Region,1))/(NoOfWindows+1));
			ww = 0;
			for ll = 0:OverlapSize:(NoOfWindows-1)/2
				bb = Speech_Region(floor(ll*WindowSize)+1:floor(ll*WindowSize)+WindowSize).*hamming(WindowSize);
				fb = fft(bb);                
				mb = 2595 * log10(1 + fb./700);                
				mfout = dct(log(abs(mb)),NoOfFilters);                    
				mfcc_test(1,ww*NoOfFilters+1:ww*NoOfFilters+NoOfFilters) = mfout;                                                        
				ww = ww + 1;
			end            
					   
			%% Classify MFCC test data on Mahanalobis distance
			D1(1) = mahal(obj_A,mfcc_test);
			D1(2) = mahal(obj_B,mfcc_test);
			D1(3) = mahal(obj_C,mfcc_test);
			D1(4) = mahal(obj_Five,mfcc_test);
			D1(5) = mahal(obj_Point,mfcc_test);
			D1(6) = mahal(obj_V,mfcc_test);
			[m Ind] = min(D1);
			
			if(Ind == ii)
                MFCCCorrect_Test(jj,ii) = MFCCCorrect_Test(jj,ii) + 1;
            end
			%% Extract LPC for test data
			Samples = eval(char(InputIn));
			zz = find(Samples) < max(Samples/3);%Threshold speech regions
			Samples(zz) = 0;
			zz = find(Samples);
			Speech_Region = Samples(zz);  
			ww = 1;
			lpc_test = lpc(Speech_Region,NoOfLPCFilters);
						
			%% Classify LPC test data on Mahanalobis distance
			D2(1) = mahal(obj_A2,lpc_test(2:end));
			D2(2) = mahal(obj_B2,lpc_test(2:end));
			D2(3) = mahal(obj_C2,lpc_test(2:end));
			D2(4) = mahal(obj_Five2,lpc_test(2:end));
			D2(5) = mahal(obj_Point2,lpc_test(2:end));
			D2(6) = mahal(obj_V2,lpc_test(2:end));
			[m Ind2] = min(D2);
			if(Ind2 == ii)
                LPCCorrect_Test(jj,ii) = LPCCorrect_Test(jj,ii) + 1;
            end
		end
	end
end
display('% LPC got correct is-');
sum(sum(LPCCorrect_Test))/(size(Letters,2)*size(Users,2)*NoOfTestSamples)
display('% MFCC got correct is-');
sum(sum(MFCCCorrect_Test))/(size(Letters,2)*size(Users,2)*NoOfTestSamples)