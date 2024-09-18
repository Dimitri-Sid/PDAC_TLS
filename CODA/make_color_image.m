function imcolor=make_color_image(im,cmap,w)
if ~exist('cmap','var')
    cmap=[121 248 252;... % 1 islet
    0    0    255;... % 2 duct
    80 237 80;...     % 3 blood vessel
    255  255  0;...   % 4 fat
    149 35  184;...   % 5 acinus
    255 194 245;...   % 6 connective tissue
    255 255 255;...   % 7 whitespace
    255  0  0;...     % 8 PanIN
    240 159 10];      % 9 PDAC 
end

if exist('w','var');im(im==0)=w;end
am=cmap(:,1);bm=cmap(:,2);cm=cmap(:,3);
imcolor=uint8(cat(3,am(im),bm(im),cm(im)));

