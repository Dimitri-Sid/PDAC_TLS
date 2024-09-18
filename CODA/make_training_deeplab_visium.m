%% Intro 
path(path,'base');
warning ('off','all');

% CHANGE pth to the ndpi and xml file (add backslash)
pth='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\';
%pth='\\fatherserverdw\ashleyex\Pancreas images\Fatemeh PDAC\annotations\endothelium\';
% CHANGE pth to the 10x tif files that you made
pthim='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\10x\'; % path to tif images to classify
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

%% create fixed testing TA images before running the model
% pth to testing annotations and tif images
pthtest=[pth,'testing data\']; 
pthtestim='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\testing data\10x\';
%pthtestim=pthim;

outpth=[pthtestim,'TA\'];
if ~isfolder(outpth);mkdir(outpth);end
imlist=dir([pthtestim,'*tif']);

for k=1:length(imlist)
    nm=imlist(k).name;
    im=imread([pthtestim,nm]);
    im=im(:,:,2);
    im=im<200;
    im=bwareaopen(im,20);
    imwrite(im,[outpth,nm]);
    
    disp([k length(imlist)])
end

%% Set up model
% date of model training
nm='10_23_2023';

% resolution of images used for training
umpix=1; % um/pixel of images used % 1=10x, 2=5x, 4=16x

% Previous attempt:
%       islet pdac normal duct ecm fiducials TLS acinar nerve fat blank space vasculature
%       islet pdac duct ecm fiducials TLS acini nerve fat white vasculature
%    #: 1     2    3    4   5         6   7     8     9   10    11         
% new#: 1     2    2    3   4         5   6     7     8   9     10
%   ws: dl    dl   dl   dl  dl        dl  dl    dl    kp  kp    dl
% nest: 
% % define actions to take per original annotation class 
% WS{1}=[0 0 0 0 0 0 0 0 2 2 0];           % remove whitespace if 0, keep only whitespace if 1, keep both if 2
% WS{2}=[10 10];                           % add removed whitespace to this class
% WS{3}=[1 2 2 3 4 5 6 7 8 9 10];         % rename classes accoring to this order 
% WS{4}=[4 6 1 2 3 5 7 8 9 11 10];       % reverse priority of classes
% WS{5}=[];                                % delete classes
% numclass=max(WS{3});
% sxy=1000;
% pthDL=[pth,nm,'\'];
% nwhite=WS{3};nwhite=nwhite(WS{2});nwhite=nwhite(1);
% nblack=numclass+1;
% 
% % CHANGE THIS: color pallete for check annotaiton images and classification
% cmap=[149 35  184;... % 1 islet
%     0    0    255;... % 2 pdac
%     80 237 80;...     % 3 ecm
%     %255 255 255;...   % 4 fiducial, white
%     255 0 0;...   % 4 fiducial, red
%     255 194 245;...   % 5 TLS
%     240 100 10;...    % 6 acini 
%     173 216 230;...     % 7 nerve 
%     240 250 10;...    % 8 fat 
%     110 90 40;...     % 9 whitespace 
%     000 000 000];     % 10 vasculature
% 
% classNames = ["islet" "pdac" "ecm" "fiducial" "TLS" "acini" "nerve" "fat" "whitespace" "vasculature" "black"];
% cmap2=cat(1,[0 0 0],cmap)/255;
% 
% titles=["islet" "pdac" "ecm" "fiducial" "TLS" "acini" "nerve" "fat" "whitespace" "vasculature"];
% make_cmap_legend(cmap,titles)

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
cmap2=cat(1,[0 0 0],cmap)/255;
titles=["islet" "PDAC" "duct" "ECM" "fiducials" "TLS" "acini" "nerve" "fat" "whitespace","vasculature"];
make_cmap_legend(cmap,titles)

ntrain=9;
nvalidate=3;
% ntrain and nvalidate in the build_model_tiles correspond to the number of “big tiles” that we are going to create for the training and the validation of the model respectively. 
% Typically I would say a good value for them would be ntrain = 7 and nvalidate = 2

%% 2 load and format annotations for each image
[ctlist,numann0]=load_xml_loop(pth,pthim,WS,umpix,nm,numclass,cmap2);
%[ctlist,numann0]=load_xml_data_loop(pth,pthim,WS,umpix,nm,numclass,cmap2);

%% make training tiles for model
build_model_tiles(pthDL,classNames,nblack,sxy,numann0,ctlist, ntrain, nvalidate)
%build_training_tiles(pthDL,classNames,nblack,sxy,numann0,ctlist)

%% build model
train_deeplab(pthDL,1:numclass+1,sxy,classNames);

%% test model (will print out the anns on test image)
pthtest='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\testing data\';
pthtestim=[pthtest,'10x\'];
load_xml_loop(pthtest,pthtestim,WS,umpix,nm,numclass,cmap2);
%load_xml_data_loop(pthtest,pthtestim,WS,umpix,nm,numclass,cmap2);

%% make confusion matrix using testing data (output perfomance results table with test data monkey)
pthtestdata=[pthtest,'data\'];
deeplab_classification_AK(pthtestim,pthDL,sxy,nm,cmap,nblack,nwhite);
pthclassifytest=[pthtestim,'classification_',nm,'\'];

%% test model performance and plot confusion matrix 
%test_data_monkey(pthtestdata,pthclassifytest,nwhite,nblack,classNames); 
pthtestdata='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\testing data\data\';
pthclassifytest='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\testing data\10x\classification_10_23_2023\';
test_model_performance(pthtestdata,pthclassifytest,nwhite,nblack,classNames(1:end-1)); 

%% classify using trained model (classify original images and new images)
% classify
pthclassify='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\10x\';
%pthclassify='\\fatherserverdw\ashleyex\Pancreas images\Fatemeh PDAC\PDAC S1\10x\';
deeplab_classification_AK(pthclassify,pthDL,sxy,nm,cmap,nblack,nwhite);

%% sanity checks
im= imread('C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\10x\classification_09_12_2023\1568s2.tif'); 
figure;imagesc(im)

