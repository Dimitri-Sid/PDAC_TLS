%% colors for classification - CHANGE to match Dimitri's classes and make legend
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


%% visualize what the volcell looks like:
%with im being the original 10x HE image

% im = imread('C:\Users\Lucie\Downloads\Dimitri\1568s1.tif');   % here it should be the path to the 10x image in your own folder
% volcell = load('C:\Users\Lucie\Downloads\cells1568s1.mat');   % same but with the path to the volcell image

im= imread('C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\cell detection\makevolcell_1568s1\1568s1.tif'); 
load('C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\cell detection\makevolcell_1568s1\cells1568s1.mat','volcell'); 
im0 = imresize(im,0.25);
imshowpair(im0,volcell);

%% set up 
% % path to classified images
% pthim='\\motherserverdw\AshleySync\other collaborations\Elana Fertig\visium paper\10x\classification_10_11_2022\'; % CHANGE
% % path to "volcell" style images (binary image with cell coordinates
% pthcell='\\motherserverdw\AshleySync\other collaborations\Elana Fertig\visium paper\10x\fix stain\Hchannel\cell_coords\cell_images\'; % CHANGE
% % place to save cropped classified images and cell coordinates
% outpth='\\motherserverdw\AshleySync\other collaborations\Elana Fertig\visium paper\10x\cropped images\'; % CHANGE
% mkdir(outpth);
% % name of image to crop
% nm='CD05B.tif'; % change this and call the code for each of the 3 images you have
% % name of output files
% outnm='CP05B'; % beginning of the filename
% % number of tissues per file
% numtiss=4;

% path to classified images
pthim='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\CODA_annotations_transfer_to_Visium\classified_images\';
% path to "volcell" style images (binary image with cell coordinates
pthcell='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\CODA_annotations_transfer_to_Visium\volcell_images\'; % DUMP THE cells*NAME*.mat FILES FROM THE makevolcell folders after renaming them to match exactly the classified images folder
% place to save cropped classified images and cell coordinates
outpth='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\CODA_annotations_transfer_to_Visium\cropped images\'; % CHANGE
mkdir(outpth);
% name of image to crop
nm='1568s1.tif'; % change this and call the code for each of the 3 images you have
% name of output files
outnm='1568s1'; % beginning of the filename
% number of tissues per file
numtiss=4;

im=imread([pthim,nm]);
imcolor=make_color_image(im,cmap,7);
% imcell=imread([pthcell,nm]);
load([pthcell,strrep(nm,'tif','mat')],'volcell');
imcell=volcell;

% resize volcell to same size as classified image, do not need if we save volcell as 10x
ii=find(volcell>0);
[x,y]=ind2sub(size(volcell),ii); 
imcell=zeros(size(im));
x=x*4;y=y*4;
ii=sub2ind(size(imcell),x,y);
imcell(ii)=1;
%figure;imshowpair(im,imcell) %check match in case xy are swapped
tmp=imgaussfilt(imcell,3);figure;imshowpair(im,tmp) %check match in case xy are swapped, this is smoother

%% crop and save, DO NOT INCLUDE TISSE OUTSIDE OF FIDUCIALS. double click to crop and move on
for z=1:numtiss
    h=figure;[~,rr]=imcrop(imcolor);close(h);
    rr=round(rr);
    
    imout=imcrop(im,rr);
    imcolorout=imcrop(imcolor,rr);
    imcellout=imcrop(imcell,rr);
    
    outpth2=[outpth,outnm,'_',num2str(z),'\'];mkdir(outpth2);
    imwrite(imout,[outpth2,'classified.tif']);
    imwrite(imcolorout,[outpth2,'color classified.tif']);
    imwrite(imcellout,[outpth2,'cell_count.tif']);
    save([outpth2,'cropvals.mat'],'rr');
end


