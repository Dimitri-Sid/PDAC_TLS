%Two subfolders, one containing Visium spaceranger data and one containing CODA outputs. 
% the name of these folders indicates the TMA being analyzed, the name of all files should be identical between TMAs
% pthS0='\\motherserverdw\AshleySync\othercollaborations\Elana Fertig\visium paper\spaceranger\';
% pth0='\\motherserverdw\AshleySync\othercollaborations\Elana Fertig\visium paper\10x\cropped images\';
% nms=dir(pth0);nms=nms(3:end);
pthS0='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\J1568_Visium\spaceRanger_output\EJ01JHU518_000_analysis\spaceranger\count\';
pth0='C:\Users\sidir\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA\CODA_annotations_transfer_to_Visium\cropped images\renamed to match spaceranger\';
nms=dir(pth0);nms=nms(3:end);


% if you classify with 10x in CODA, resolution = 1 micron
% if you classify with 5x in CODA, resolution = 2 micron

%%
nns=ones([1 12])*0.26; % multiple visium coords to scale and match the cropped images  MATCHES THE XCEL DATA TO THE VISIUM IMAGE
%nns(6:12)=0.26; %& need to change values accordingly so each image has the right scaling factor, need to visualize manually 
for k=1:length(nms)
    pth=[pth0,nms(k).name,'\'];
    pthS=[pthS0,nms(k).name,'\spatial\'];
    %pthS=[pthS0,'spaceranger_output_',nms(k).name,'\'];
    im=imread([pthS,'tissue_hires_image.png']);
    imC=double(imread([pth,'classified.tif']));
    imCcolor=imread([pth,'color classified.tif']);
    imcell=double(imread([pth,'cell_count.tif'])>0);
    imF=imC==5;imF=bwareaopen(imF,500); % fiducial markers
    
    % rename classification to combine whitespace with fiducial, base on the new#
        %my annotation scheme:
    % %       islet pdac normal duct ecm fiducials TLS acinar nerve fat blank space vasculature
    % %       islet pdac duct ecm fiducials TLS acini nerve fat white vasculature
    % %    #: 1     2    3    4   5         6   7     8     9   10    11        

    imC(imC==10)=5; %imC is the classified image, showing it in new#
    imC(imC==11)=10;
    imcell(imC==5)=0; %imcell has the nuclei, whatever is the whitespace/fiducials now you dont count the nuclei in it to remove fasle positives
    imcell=imcell.*imC; %masking operation, the.* tells matlab to multiply 2 matrices 
    
    nmX='tissue_positions.csv';
    [~,~,data]=xlsread([pthS,nmX]);
    outfile=[pth,nmX(1:end-4),''];
    nn=nns(k);
    titles={'','','','','','','islet','pdac','duct','ecm','whitespace/fiducials','TLS', 'acinar','nerve','fat','vasculature'}; % this is the updates titles that corresponds to the new classes we made(new2#), need to keep the first 6 empty

    im = im2uint8(im); 
    data = data(2:end,:);
    register_visium_to_classified_image_cell_count(im,imC,imCcolor,imcell,data,outfile,imF,nn,titles);
    close all
    disp([k length(nms)]);
end




