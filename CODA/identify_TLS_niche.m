
% 4 means it is closing within 4 pixels 
% ANYTHING WHITE KEEPING
% ANYTHING GREEN REMOVING
% RED IS PIXELS GAINED WITH CLOSE

%pth='C:\Users\akiemen1\Documents\'; % THIS SHOULD BE THE PATH TO THE 10x CLASSSIFIED IMAGES
%pthHE='C:\Users\akiemen1\Documents\HE\'; % THIS SHOULD BE THE PATH TO THE 10x H&E IMAGES
%pthCC='\'; % THIS SHOULD BE THE PATH TO THE CHECK ANNOTATION IMAGES

pth='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\J1568_HE_Slides_Procured_8_2022\10-25-22 j1568 annotation edits\10x\classification_10_23_2023\'; % THIS SHOULD BE THE PATH TO THE 10x CLASSSIFIED IMAGES
pthHE='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\J1568_HE_Slides_Procured_8_2022\10-25-22 j1568 annotation edits\10x\'; % THIS SHOULD BE THE PATH TO THE 10x H&E IMAGES
pthCC='C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\J1568_HE_Slides_Procured_8_2022\10-25-22 j1568 annotation edits\10x\classification_10_23_2023\check_classification\';
opath = 'C:\Users\dsidiro1\OneDrive - Johns Hopkins\FertigLab\Spatial\J1568_trial\CODA whole slide\metanalyses\';

imlist=dir([pth,'*tif']);

minTLS=5000; % delete TLS smaller than this many pixels - YOU CAN TUNE THIS
maxTLS=50000; % delete TLS larger than this many pixels - YOU CAN TUNE THIS
size_window=300; % pixel size in 10x - YOU CAN TUNE THIS

downsample=2; % downsample images to speed up calculation
size_window=round(size_window/downsample);

TLS=6;
WS=10;

titles={"islet","PDAC","duct","ecm","fiducials","TLS","acini","nerve","fat","nontissue","vasculature"};

for k=1:length(imlist)
    nm=imlist(k).name;
    im0=imread([pth,nm]);
    im=im0(1:downsample:end,1:downsample:end,:); % quantify at 5x to speed up calculations
    outpth=[opath,'TLS_data_',num2str(size_window), '_',num2str(minTLS),'_',num2str(maxTLS),'\',strrep(nm,'.tif','\')];mkdir(outpth);

    % extract TLS and smooth
    imTLS0=im==TLS;
    imTLS=imclose(imTLS0,strel('disk',7));
    imTLS=imfill(imTLS,'holes');
    imTLS=bwareaopen(imTLS,minTLS);        % delete pixels less than XX
    imTLS=imTLS-bwareaopen(imTLS,maxTLS); % delete pixels more than XX
    figure(15);imshowpair(imTLS0,imTLS)  % white = things we kept, green = things we deleted bc too small or too large
    
    % for each TLS, extract it's neighboorhood
    TLSmask=bwlabel(imTLS); % label each TLS independently
    TLScomp=zeros([max(TLSmask(:)) 11]);
    for b=1:max(TLSmask(:))
        tmp=TLSmask==b;
        im_window0=bwdist(tmp)<=size_window;
        im_window=im_window0.*~tmp;
        im_window=im_window.*double(im);
        comp=histcounts(im_window(:),1:12);
        comp(WS)=0;comp=comp/sum(comp)*100;
        TLScomp(b,:)=comp;

        a1=sum(im_window0,2);a2=sum(im_window0,1);
        bbHE=[find(a1>0,1,'first') find(a1>0,1,'last') find(a2>0,1,'first') find(a2>0,1,'last')]*2;
        HEtile=imread([pthHE,nm],'PixelRegion',{[bbHE(1) bbHE(2)],[bbHE(3) bbHE(4)]});
        
        % save H&E image
        outnm=strrep(nm,'.tif',['_TLS_',num2str(b),'.tif']);
        imwrite(HEtile,[outpth,outnm])
        
        % save check annotations image
        bb=bbHE/2;
        imcheck=imread([pthCC,strrep(nm,'.tif','.jpg')]);
        imcheck=imresize(imcheck,size(im),'nearest');
        imcheck=imcheck(bb(1):bb(2),bb(3):bb(4),:);
        imcheck=imresize(imcheck,size(HEtile(:,:,1)),'nearest');
        outnmCC=strrep(outnm,'.tif','_check_annotation.tif');
        imwrite(imcheck,[outpth,outnmCC]);

    end
    filename=[outpth,strrep(nm,'.tif','.xlsx')];
    T=cat(1,titles,num2cell(TLScomp));
    writecell(T,filename,'Sheet',1)

end