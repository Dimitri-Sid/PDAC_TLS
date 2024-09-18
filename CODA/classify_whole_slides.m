%% Intro 
path(path,'base');
warning ('off','all');

% CHANGE pth to the ndpi and xml file (add backslash)
pth='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA whole slide\';
%pth='\\fatherserverdw\ashleyex\Pancreas images\Fatemeh PDAC\annotations\endothelium\';
% CHANGE pth to the 10x tif files that you made
pthim='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\J1568_HE_Slides_Procured_8_2022\10-25-22 j1568 annotation edits\10x\'; % path to tif images to classify
%pthim='\\fatherserverdw\ashleyex\Pancreas images\Fatemeh PDAC\annotations\10x\'; % path to tif images to classify

%% create fixed TA images from training images before running the model
% from fix_TA_images.m
% ashley on 10/5: call the attached code “fix_TA_images” to make bett er TA files that more accurately detect the ti ssue space in the 10x H&E images.
% Below, the top row is the original TA image and the bottom row is using the improved code. The goal of this image is to detect ti ssue space so that we can automati cally remove whitespace from the pdac, vasculature, ecm, etc..annotati ons. When the TA was messed up we kept whitespace in those annotati ons, resulti ng in the misclassifi cati on of whitespace as ti ssue that you noted in your last model. This should fi x that problem
% pth='other collaborations\Elana Fertig\Dimitri\imagescope files\10x\'; % REPLACE THIS WITH THE PATH TO YOUR 10X IMAGES
outpth=[pthim,'TA\'];
if ~isfolder(outpth);mkdir(outpth);end
imlist=dir([pthim,'*tif']);

for k=1:length(imlist)
    nm=imlist(k).name;
    im=imread([pthim,nm]);
    im=im(:,:,2);
    im=im<200;
    im=bwareaopen(im,20);
    imwrite(im,[outpth,nm]);
    
    disp([k length(imlist)])
end


%% Seelect model
% date of model training
nm='10_23_2023';

% resolution of images used for training
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% New cmap and WS code:
% define actions to take per annotation class
%  islet pdac duct ecm fiducials TLS acini nerve fat white vasculature
WS{1}=[2 0 0 0 0 2 2 2 2 2 0];    % remove whitespace if 0, keep only whitespace if 1, keep both if 2
WS{2}=[10 4];                     % add removed whitespace to this class
WS{3}=[1 2 3 4 5 6 7 8 9 10 11];  % rename classes accoring to this order
WS{4}=[4 7 6 9 1 2 3 8 11 10 5];  % reverse priority of classes
WS{5}=[];                         % delete classes
numclass=max(WS{3});
sxy=1000;
pthDL=[pth,nm,'\'];
nblack=numclass+1;
nwhite=WS{2};nwhite=nwhite(1);

%% make legend
cmap=[121 248 252;... % 1  islet
    240 159 010;...   % 2  PDAC
    000 000 255;...   % 3  duct
    255 194 245;...   % 4  ecm
    125 125 125;...   % 5  fiducials
    000 000 000;...   % 6  TLS
    149 035 184;...   % 7  acini
    073 120 111;...   % 8  nerve
    255 255 000;...   % 9  fat
    255 255 255;...   % 10 whitespace
    080 237 080];     % 11 blood vessel / vasculature

classNames = ["islet" "PDAC" "duct" "ECM" "fiducials" "TLS" "acini" "nerve" "fat" "whitespace" "vasculature" "black"];
titles=["islet" "PDAC" "duct" "ECM" "fiducials" "TLS" "acini" "nerve" "fat" "whitespace","vasculature"];
make_cmap_legend(cmap,titles)

%% classify using trained model (classify original images and new images)
% classify
pthclassify=pthim;
%pthclassify='\\fatherserverdw\ashleyex\Pancreas images\Fatemeh PDAC\PDAC S1\10x\';
deeplab_classification_AK(pthclassify,pthDL,sxy,nm,cmap,nblack,nwhite);

%% sanity checks
im= imread('C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\10x\classification_09_12_2023\1568s2.tif'); 
figure;imagesc(im)

