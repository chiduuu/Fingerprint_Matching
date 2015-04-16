function GUI

%SPLASH(FILENAME,FMT,TIME)
 SPLASH('Fprint_logo','bmp',3000);
 
hMainFigure = figure(...     % The main GUI figure
                    'Name','Image Tranforms','NumberTitle','off',...
                    'MenuBar','none', ...
                    'Toolbar','none', ...
                    'Color', get(0,...
                             'defaultuicontrolbackgroundcolor'));                         
hFileMenu      =   uimenu(...       % File menu
                        'Parent',hMainFigure,...
                        'HandleVisibility','callback', ...
                        'Label','File');
                   
hOpenMenuitem  =   uimenu(...       % Open menu item
                        'Parent',hFileMenu,...
                        'Label','Open',...
                        'HandleVisibility','callback', ...
                        'Callback', @hOpenMenuitemCallback);  
                    
hCleareMenuitem  =  uimenu(...       % Close menu item
                        'Parent',hFileMenu,...
                        'Label','Clear',...
                        'Separator','on',...
                        'HandleVisibility','callback', ...
                        'Callback', @hClearMenuitemCallback);                    
                  
hCloseMenuitem  =  uimenu(...       % Close menu item
                        'Parent',hFileMenu,...
                        'Label','Close',...                        
                        'HandleVisibility','callback', ...
                        'Callback', @hCloseMenuitemCallback);
                    
hFunctionMenu      =   uimenu(...       % Function menu
                        'Parent',hMainFigure,...
                        'HandleVisibility','callback', ...
                        'Label','Function');
                    
hCorePointMenuitem  =   uimenu(...       % DWT menu item
                        'Parent',hFunctionMenu,...
                        'Label','Core Point',...
                        'HandleVisibility','callback', ...
                        'Enable','off', ...
                        'Callback', @hcorepointCallback);
                    
hCropMenuitem  =   uimenu(...       % core point menu item
                        'Parent',hFunctionMenu,...
                        'Label','Crop Image',...
                        'HandleVisibility','callback', ...
                        'Enable','off', ...
                        'Callback', @hcropCallback); 
                    
hTransformMenu      =   uimenu(...       % Transform menu
                        'Parent',hMainFigure,...
                        'HandleVisibility','callback', ...
                        'Label','Transform');
                    
hDWTMenuitem  =   uimenu(...       % DWT menu item
                        'Parent',hTransformMenu,...
                        'Label','DWT',...
                        'Enable','off', ...
                        'HandleVisibility','callback'); 
                    
hDWTl1Menuitem  =   uimenu(...       % DWT menu item
                        'Parent',hDWTMenuitem,...
                        'Label','DWT Level 1',...
                        'HandleVisibility','callback', ...
                        'Callback', @hDWTl1MenuitemCallback); 
hDWTl2Menuitem  =   uimenu(...       % DWT menu item
                        'Parent',hDWTMenuitem,...
                        'Label','DWT Level 2',...
                        'HandleVisibility','callback', ...
                        'Callback', @hDWTl2MenuitemCallback); 
hDWTl3Menuitem  =   uimenu(...       % DWT menu item
                        'Parent',hDWTMenuitem,...
                        'Label','DWT Level 3',...
                        'HandleVisibility','callback', ...
                        'Callback', @hDWTl3MenuitemCallback);
hDWTFVMenuitem  =   uimenu(...       % DWT menu item
                        'Parent',hDWTMenuitem,...
                        'Label','DWT Feature vector',...
                        'HandleVisibility','callback', ...
                        'Callback', @hDWTFVMenuitemCallback);                     

hDatabaseMenu      =   uimenu(...       % Database menu
                        'Parent',hMainFigure,...
                        'HandleVisibility','callback', ...
                        'Label','Database'); 

hAdditem  =   uimenu(...       %  add item
                        'Parent',hDatabaseMenu,...
                        'Label','Add Image',...
                        'HandleVisibility','callback', ...
                        'Callback', @hAdditemCallback);       
                    
hDelitem  =   uimenu(...       % Del menu item
                        'Parent',hDatabaseMenu,...
                        'Label','Delete Database',...
                        'HandleVisibility','callback', ...
                        'Callback', @hDelitemCallback); 
                    
hMatchMenuitem  =   uimenu(...       % Match menu item
                        'Parent',hFunctionMenu,...
                        'Label','Match',...
                        'HandleVisibility','callback',...
                        'Enable','off', ...
                        'Callback', @hMatchitemCallback); 

hHelpMenu      =   uimenu(...       % Help menu
                        'Parent',hMainFigure,...
                        'HandleVisibility','callback', ...
                        'Label','Help');
                    
hAboutitem  =   uimenu(...       % Match menu item
                        'Parent',hHelpMenu,...
                        'Label','About us',...
                        'HandleVisibility','callback',...
                        'Callback', @hAboutitemCallback);   
          
hContactitem  =   uimenu(...       % Contact menu item
                        'Parent',hHelpMenu,...
                        'Label','Contact us',...
                        'HandleVisibility','callback',...
                        'Callback', @hContactitemCallback); 
                    
%%%%%%%%%%%%%%%%%%%%%Core Point function%%%%%%%%%%%%%%%%%%%%
function [XofCenter,YofCenter] = centralizing(fingerprint)
global immagine n_bands h_bands n_arcs h_radius h_lato n_sectors matrice

x=[-16:1:16];
y=[-16:1:16];
%size of the dimension of X & Y specified by scalar dim
dimx=size(x,2);
dimy=size(y,2);
% variance gaussiana, order filter complex
variance=sqrt(55);
order=1;
gamma=2;
filter_core=zeros(dimx,dimy);
filter_delta=zeros(dimx,dimy);
for ii=1:dimx
    for jj=1:dimy
        exponent=exp(-(x(ii)^2+y(jj)^2)/(2*variance^2));
        % filter core
        filter_factor=x(ii)+i*y(jj);
        filter_core(ii,jj)=exponent*filter_factor^order;
        % filter delta        
        filter_factor=x(ii)-i*y(jj);
        filter_delta(ii,jj)=exponent*filter_factor^order;
    end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------------------------------------
% The low-pass filter ---------------
%------------------------------------
% Gaussian Low Pass Filter ----------
%------------------------------------
x=[-16:1:16];
y=[-16:1:16];
dimx=size(x,2);
dimy=size(y,2);
variance=sqrt(1.2);
filter=zeros(dimx,dimy);
for ii=1:dimx
    for jj=1:dimy
        exponent=exp(-(x(ii)^2+y(jj)^2)/(2*variance^2));
        filter(ii,jj)=exponent;
    end
end
% normalization
filter=filter/sum(sum(filter)); 
%------------------------------------
%------------------------------------
img=fingerprint;
img=double(img);
%--------------------------------------------------------------------------
% complex field at 0 level
[gx,gy]=gradient(img);
num=(gx+i*gy).^2;
den=abs((gx+i*gy).^2);
pos=find(den);
num(pos)=num(pos)./den(pos);
z=zeros(size(img,1),size(img,2));
z(pos)=num(pos);
pos=find(den==0);
z(pos)=1;
%**********************************
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if 1==0%------pyramid gaussian
    %-------------------------------------
    % complex(orientation) field at  level 1 and call to conv2fft file
    z1=conv2fft(z,filter,'same');
    z1=dyaddown(z1,1,'m');
    num=z1;
    den=abs(z1);
    pos=find(den);
    num(pos)=num(pos)./den(pos);
    z1=zeros(size(z1,1),size(z1,2));
    z1(pos)=num(pos);
    pos=find(den==0);
    z1(pos)=1;
    %**********************************
    %-------------------------------------
    % complex field at  level 2 and call to conv2fft file
    z2=conv2fft(z1,filter,'same');
    z2=dyaddown(z2,1,'m');
    num=z2;
    den=abs(z2);
    pos=find(den);
    num(pos)=num(pos)./den(pos);
    z2=zeros(size(z2,1),size(z2,2));
    z2(pos)=num(pos);
    pos=find(den==0);
    z2(pos)=1;
    %**********************************
    %-------------------------------------
    % complex field at  level 3 and call to conv2fft file
    z3=conv2fft(z2,filter,'same');
    z3=dyaddown(z3,1,'m');
    num=z3;
    den=abs(z3);
    pos=find(den);
    num(pos)=num(pos)./den(pos);
    z3=zeros(size(z3,1),size(z3,2));
    z3(pos)=num(pos);
    pos=find(den==0);
    z3(pos)=1;
    %**********************************
    %-------------------------------------
    %-----------------------------------------> z z1 z2 z3 -----------------------
    % complex field per each levels -----------------------
    %-------------------------------------------------------------------
    z_f=conv2fft(z,filter_core,'same');
    z_f=abs(z_f);
    temp0=z_f;%temp0
    %-------------------
    z_1f=conv2fft(z1,filter_core,'same');
    z_1f=abs(z_1f);
    temp1=z_1f;%temp1
    %-------------------
    z_2f=conv2fft(z2,filter_core,'same');
    z_2f=abs(z_2f);
    temp2=z_2f;%temp2
    %---------------------
    z_3f=conv2fft(z3,filter_core,'same');
    z_3f=abs(z_3f);
    temp3=z_3f;%temp3
    %-------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %------------------------------------
    % detailed ----------------
    %------------------------------------
    % usage of temp0 temp1 temp2 temp3
    %------------------------------------
    [maximus_vector,position_vector]=max(temp3);
    [maximus,position]=max(maximus_vector);
    y_max=position;
    x_max=position_vector(position);
    maximus;
    
    x0=2*x_max;
    y0=2*y_max;
    
    dx=10;
    dy=10;
    
    positioni=zeros(size(temp2));
    positioni(max(1,x0-dx):min(size(temp2,1),x0+dx),max(1,y0-dy):min(size(temp2,2),y0+dy))=1;
    temp2=temp2.*positioni;
    
    [maximus_vector,position_vector]=max(temp2);
    [maximus,position]=max(maximus_vector);
    y_max=position;
    x_max=position_vector(position);
    maximus;
    
    x0=2*x_max;
    y0=2*y_max;
    
    dx=10;
    dy=10;
    
    positioni=zeros(size(temp1));
    positioni(max(1,x0-dx):min(size(temp1,1),x0+dx),max(1,y0-dy):min(size(temp1,2),y0+dy))=1;
    temp1=temp1.*positioni;
    
    [maximus_vector,position_vector]=max(temp1);
    [maximus,position]=max(maximus_vector);
    y_max=position;
    x_max=position_vector(position);
    maximus;
    
    x0=2*x_max;
    y0=2*y_max;
    
    dx=5;
    dy=5;
    
    positioni=zeros(size(temp0));
    positioni(max(1,x0-dx):min(size(temp0,1),x0+dx),max(1,y0-dy):min(size(temp0,2),y0+dy))=1;
    temp0=temp0.*positioni;
    
    [maximus_vector,position_vector]=max(temp0);
    [maximus,position]=max(maximus_vector);
    y_max=position;
    x_max=position_vector(position);
    maximus;
    
    disp('Coordinate x y');
    disp(x_max);
    disp(y_max);
    
    XofCenter=y_max;
    YofCenter=x_max;
    Outputprint=zeros(50);
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%-------------------------------- -------------------------------------->
%                          ---------- -----------------------  ------------ -             
%-----------------------------------------> z z1 z2 z3 -----------------------
% complex field for different levels -----------------------
%-------------------------------------------------------------------
%--------------------------------------------------------------------------
% ------------------- parameters -------------------------------------------
angle=0;        % angle rotationise immagine initialise
bxv=8;          % dimensions block variance
byv=8;
bxc=64;         % dimensions block core point
byc=64;
threshold_var=20;  % threshold variance
dimseclose=10;  % dimension operation closing
dimseerode=44;  % dimension  operation erosion
maxcore=200;   % max number of core points calculated during scanning
[dimx,dimy]=size(fingerprint);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
temp=z;
[temp,dimxt,dimyt]=mirror(temp);
z_f=conv2fft(temp,filter_core,'same');
z_f=recrop(z_f,dimxt,dimyt);
z_f=abs(z_f);
%--------------------------------------------
%--------------------------------------------
% resize-------------------- 
imgd=double(fingerprint);
dimxr=dimx-mod(dimx,bxv);
dimyr=dimy-mod(dimy,byv);
imgr=imgd(1:dimxr,1:dimyr);
%---------------------------
nbx=dimxr/bxv;
nby=dimyr/byv;
mat_var=zeros(dimxr,dimyr);
for ii=1:nbx
    for jj=1:nby
        block=imgr((ii-1)*bxv+1:ii*bxv,(jj-1)*byv+1:jj*byv);
        media=sum(sum(block))/(bxv*byv);
        variance=1/(bxv*byv)*sum(sum(abs(media.^2-block.^2)));
        mat_var((ii-1)*bxv+1:ii*bxv,(jj-1)*byv+1:jj*byv)=sqrt(variance);
    end
end
mat_ok=zeros(dimxr,dimyr);
pos=find(mat_var>threshold_var);
mat_ok(pos)=1;
mat_ok(dimx,dimy)=0;
%figure('Name','variance simplified > threshold');
%imshow(mat_ok);

mat_ok=imclose(mat_ok,ones(dimseclose));
%figure('Name','variance Closed');
%imshow(mat_ok);

mat_ok=imerode(mat_ok,ones(dimseerode));
%figure('Name','variance Closed ed Eroded');
%imshow(mat_ok);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% calculate core point at every block
dimxr=dimx-mod(dimx,bxc);
dimyr=dimy-mod(dimy,byc);
imgr=imgd(1:dimxr,1:dimyr);
matrice_finale=z_f.*mat_ok;
%--------------------------------------------------------------------------
%------------------------------ not used-----------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if 1==0%--------------------------------> un maximus per the block
    nbx=dimxr/bxc;
    nby=dimyr/byc;
    mat_core=zeros(maxcore,2);
    count=1;
    for ii=1:nbx
        for jj=1:nby
            block = matrice((ii-1)*bxc+1:ii*bxc,(jj-1)*byc+1:jj*byc);
            [maximus_vector,position_vector]=max(block);
            [maximus,position]=max(maximus_vector);
            y_max=position;
            x_max=position_vector(position);
            x0=(ii-1)*bxc+1+(x_max-1);
            y0=(jj-1)*byc+1+(y_max-1);
            if maximus>0
                mat_core(count,1)=x0;
                mat_core(count,2)=y0;
                count=count+1;
            end
        end
    end
    out=mat_core;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

[maximus_vector,position_vector]=max(matrice_finale);
[maximus,position]=max(maximus_vector);
y_max=position;
x_max=position_vector(position);

XofCenter=y_max;
YofCenter=x_max;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hOpenMenuitemCallback(hObject, eventdata, handles)
% Callback function run when the Open menu item is selected
% hObject    handle to mopen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.bmp','All Image Files' });
PathName=[PathName,FileName];
%msgbox(PathName)
if ~isequal(FileName, 0) 
    
clf(hMainFigure);
set(hCorePointMenuitem,'Enable','off');
set(hMatchMenuitem,'Enable','off');
set(hCropMenuitem,'Enable','off');
set(hDWTMenuitem,'Enable','off');
    
hlabel1 = uicontrol(hMainFigure,'Style','text',...
                'String','Input Image',...
                'Units', 'normalized', ...
                'Position',[0.15 0.8 0.1 0.1]);    

himageaxes = axes(...    % Axes for displaying the image
                 'Parent', hMainFigure, ...
                 'Units', 'normalized', ...
                 'HandleVisibility','callback', ...
                 'Position',[0.0 0.4 0.40 0.40]);
set(hMainFigure,'CurrentAxes',himageaxes)
img=imread(PathName);
imshow(img);
setappdata(hOpenMenuitem,'image',img); 
set(hCorePointMenuitem,'Enable','on');
set(hMatchMenuitem,'Enable','on');
end 

end

function hClearMenuitemCallback(hObject, eventdata, handles)
   selection = ...
       questdlg(['Clear ' get(hMainFigure,'Name') '?'],...
                ['Clear ' get(hMainFigure,'Name') '...'],...
                'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end 
    clf(hMainFigure);
    set(hCorePointMenuitem,'Enable','off');
    set(hMatchMenuitem,'Enable','off');
    set(hCropMenuitem,'Enable','off');
    set(hDWTMenuitem,'Enable','off');
end

function hCloseMenuitemCallback(hObject, eventdata, handles)
% Callback function run when the Close menu item is selected
    selection = ...
       questdlg(['Close ' get(hMainFigure,'Name') '?'],...
                ['Close ' get(hMainFigure,'Name') '...'],...
                'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end

    delete(hMainFigure);
end

function hcorepointCallback(hObject, eventdata, handles)
% Callback function run to detect the core point

fingerprint=getappdata(hOpenMenuitem,'image');
%fingerprint=imread('101_1.tif');
[XofCenter,YofCenter] = centralizing(fingerprint);

setappdata(hCorePointMenuitem ,'i',XofCenter);
setappdata(hCorePointMenuitem ,'j',YofCenter);

hlabel2= uicontrol(hMainFigure,'Style','text',...
                'String','Core Point',...
                'Units', 'normalized', ...
                'Position',[0.45 0.8 0.1 0.1]);
        
himageaxes = axes(...    % Axes for displaying the image
                 'Parent', hMainFigure, ...
                 'Units', 'normalized', ...
                 'HandleVisibility','callback', ...
                 'Position',[0.3 0.4 0.40 0.40]); %position [left bottom width height]
set(hMainFigure,'CurrentAxes',himageaxes)
imshow(fingerprint);
hold on;
plot(XofCenter,YofCenter,'o');     %410 block
%plot(224,305,'o');
hold off;
set(hCropMenuitem,'Enable','on');
end

function hcropCallback(hObject, eventdata, handles)

img=getappdata(hOpenMenuitem,'image');
i=getappdata(hCorePointMenuitem ,'i');
j=getappdata(hCorePointMenuitem ,'j');
i=i-31;
j=j-31;
ci=imcrop(img,[i j 63 63]);
imwrite(ci,'image.tif');
%whos ci

hlabel3= uicontrol(hMainFigure,'Style','text',...
                'String','Cropped Image',...
                'Units', 'normalized', ...
                'Position',[0.75 0.8 0.1 0.1]);
        
himageaxes = axes(...    % Axes for displaying the image
                 'Parent', hMainFigure, ...
                 'Units', 'normalized', ...
                 'HandleVisibility','callback', ...
                 'Position',[0.62 0.4 0.40 0.40]); %position [left bottom width height]
set(hMainFigure,'CurrentAxes',himageaxes);
imshow(ci);
setappdata(hCropMenuitem,'image',ci); 
set(hDWTMenuitem,'Enable','on');
end

function hDWTl1MenuitemCallback(hObject, eventdata, handles)
img=getappdata(hCropMenuitem,'image');
%[cA,cH,cV,cD] = dwt2(img,'db9');
[C,S] = wavedec2(img,1,'db9');

A1 = appcoef2(C,S,'db9',1);
H1 = detcoef2('h',C,S,1);
V1 = detcoef2('v',C,S,1);
D1 = detcoef2('d',C,S,1);
%fileID = fopen('coefficients.txt','w');
%fprintf(fileID,'%d',A1,H1,V1,D1);
%whos cA
%A1=A1(1:64,1:64); 

hMainFigure1 = figure(...     % The main GUI figure
                    'Name','DWT Level 1','NumberTitle','off',...
                    'Color', get(0,...
                             'defaultuicontrolbackgroundcolor'));
                            
colormap(gray);
NBCOL = size(img,1);
tiledImage = wcodemat(A1,NBCOL);

tiledImage = [tiledImage wcodemat(H1,NBCOL); ... 
              wcodemat(V1,NBCOL) wcodemat(D1,NBCOL)]; 


image(tiledImage);
axis off          % Remove axis ticks and numbers
axis image        % Set aspect ratio to obtain square pixels

hold on
rectangle('Position',[1 1 79 79], 'LineWidth',3, 'EdgeColor','b');
line([40,40], [1,80],'LineWidth',3);  %[y1 y2], [x1 x2] 
line([1,80], [40 40],'LineWidth',3);

hold off


end


function hDWTl2MenuitemCallback(hObject, eventdata, handles)
img=getappdata(hCropMenuitem,'image');
%[cA,cH,cV,cD] = dwt2(img,'db9');

[C,S] = wavedec2(img,2,'db9');

cA2 = appcoef2(C,S,'db9',2);
cH2 = detcoef2('h',C,S,2); 
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);


NBCOL = size(img,1); 

tiledimage=[wcodemat(cA2,NBCOL) wcodemat(cH2,NBCOL); ...
            wcodemat(cV2,NBCOL) wcodemat(cD2,NBCOL)];
    
tiledimage=tiledimage(9:end-8,9:end-8);   
  
A1 = appcoef2(C,S,'db9',1);
H1 = detcoef2('h',C,S,1);
V1 = detcoef2('v',C,S,1);
D1 = detcoef2('d',C,S,1);

tiledimage=[tiledimage wcodemat(H1,NBCOL); ...
            wcodemat(V1,NBCOL) wcodemat(D1,NBCOL)];

hMainFigure1 = figure(...     % The main GUI figure
                    'Name','DWT Level 2','NumberTitle','off',...
                    'Color', get(0,...
                             'defaultuicontrolbackgroundcolor'));
                            
colormap(gray);
image(tiledimage)
axis off          % Remove axis ticks and numbers
axis image        % Set aspect ratio to obtain square pixels
 
hold on
rectangle('Position',[1 1 79 79], 'LineWidth',3, 'EdgeColor','b');
line([40,40], [1,80],'LineWidth',3);
line([1,80], [40 40],'LineWidth',3);
%%%%%
line([20,20], [1,40],'LineWidth',3);
line([1,40], [20 20],'LineWidth',3);
hold off

end

function hDWTl3MenuitemCallback(hObject, eventdata, handles)
img=getappdata(hCropMenuitem,'image');
%[cA,cH,cV,cD] = dwt2(img,'db9');

[C,S] = wavedec2(img,3,'db9');

cA3 = appcoef2(C,S,'db9',3);
cH3 = detcoef2('h',C,S,3); 
cV3 = detcoef2('v',C,S,3);
cD3 = detcoef2('d',C,S,3);
%whos cA3
NBCOL = size(img,1); 

tiledimage=[wcodemat(cA3,NBCOL) wcodemat(cH3,NBCOL); ...
            wcodemat(cV3,NBCOL) wcodemat(cD3,NBCOL)];
    
tiledimage=tiledimage(13:end-12,13:end-12);  

cA2 = appcoef2(C,S,'db9',2);
cH2 = detcoef2('h',C,S,2); 
cV2 = detcoef2('v',C,S,2);
cD2 = detcoef2('d',C,S,2);

tcH2=cH2(5:end-4,5:end-4);
tcV2=cV2(5:end-4,5:end-4);
tcD2=cD2(5:end-4,5:end-4);

NBCOL = size(img,1); 

tiledimage=[tiledimage wcodemat(tcH2,NBCOL); ...
            wcodemat(tcV2,NBCOL) wcodemat(tcD2,NBCOL)];
    
%tiledimage=tiledimage(9:end-8,9:end-8);   
  
A1 = appcoef2(C,S,'db9',1);
H1 = detcoef2('h',C,S,1);
V1 = detcoef2('v',C,S,1);
D1 = detcoef2('d',C,S,1);

tiledimage=[tiledimage wcodemat(H1,NBCOL); ...
            wcodemat(V1,NBCOL) wcodemat(D1,NBCOL)];

hMainFigure1 = figure(...     % The main GUI figure
                    'Name','DWT Level 3','NumberTitle','off',...
                    'Color', get(0,...
                             'defaultuicontrolbackgroundcolor'));
                            
colormap(gray);


image(tiledimage)
axis off          % Remove axis ticks and numbers
axis image        % Set aspect ratio to obtain square pixels


hold on
rectangle('Position',[1 1 79 79], 'LineWidth',3, 'EdgeColor','b');
line([40,40], [1,80],'LineWidth',3);
line([1,80], [40 40],'LineWidth',3);
%%%%
line([20,20], [1,40],'LineWidth',3);
line([1,40], [20 20],'LineWidth',3);
%%%%
line([10,10], [1,20],'LineWidth',3);
line([1,20], [10 10],'LineWidth',3);
hold off

end

function hDWTFVMenuitemCallback(hObject, eventdata, handles)
img=getappdata(hCropMenuitem,'image');    
fv=cell(1,36);
[cA1,cH1,cV1,cD1] = dwt2(img,'db9');
[cA2,cH2,cV2,cD2] = dwt2(cA1,'db9');
[cA3,cH3,cV3,cD3] = dwt2(cA2,'db9');
%%%extract features
%%Level 3
[fv{1} fv{2} fv{3} fv{4}]=dwt2(cH3,'db9');
[fv{5} fv{6} fv{7} fv{8}]=dwt2(cV3,'db9');
[fv{9} fv{10} fv{11} fv{12}]=dwt2(cD3,'db9');
%%Level 2
[fv{13} fv{14} fv{15} fv{16}]=dwt2(cH2,'db9');
[fv{17} fv{18} fv{19} fv{20}]=dwt2(cV2,'db9');
[fv{21} fv{22} fv{23} fv{24}]=dwt2(cD2,'db9');
%%Level 1
[fv{25} fv{26} fv{27} fv{28}]=dwt2(cH1,'db9');
[fv{29} fv{30} fv{31} fv{32}]=dwt2(cV1,'db9');
[fv{33} fv{34} fv{35} fv{36}]=dwt2(cD1,'db9');
%%%Standard Deviation
for i=1:36
      fv{i}=std(fv{i}(:));
end

%{
t=mean(fv{1}(:));
sd=0;
for i=1:19
    for j=1:19
        sd=sd+(fv{1}(i,j)-t)*(fv{1}(i,j)-t);
    end
end
fv{1}=sqrt(sd/(19*19));
%}

figure('Position', [100, 300, 600, 460],...
       'Name', 'Feature Vector',...  % Title figure
       'NumberTitle', 'off',... % Do not show figure number
       'MenuBar', 'none');      % Hide standard menu bar menus
%count = load('count.dat');
fv=fv';
tablesize = size(fv); 
% Define parameters for a uitable (col headers are fictional)
colnames = {'Standard Deviation'};
% All column contain numeric data (integers, actually)
colfmt = {'numeric'};
% Disallow editing values (but this can be changed)
coledit = [false];
% Set columns all the same width (must be in pixels)
colwdt = {150};

htable = uitable('Units', 'normalized',...
                 'Position', [0.3 0.03 0.375 0.92],...
                 'Data',  fv,... 
                 'ColumnName', colnames,...
                 'ColumnFormat', colfmt,...
                 'ColumnWidth', colwdt,...
                 'ColumnEditable', coledit);   
            


end

function hAdditemCallback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.bmp','All Image Files' });
PathName=[PathName,FileName];
if ~isequal(FileName, 0)
fingerprint=imread(PathName);
cd('Database/');
imwrite(fingerprint,FileName);
cd('../');
%copyfile(PathName,'Database/','f');
[XofCenter,YofCenter] = centralizing(fingerprint);
i=XofCenter-31;
j=YofCenter-31;
img=imcrop(fingerprint,[i j 63 63]);
fv=cell(1,36);
[cA1,cH1,cV1,cD1] = dwt2(img,'db9');
[cA2,cH2,cV2,cD2] = dwt2(cA1,'db9');
[cA3,cH3,cV3,cD3] = dwt2(cA2,'db9');
%%%extract features
%%Level 3
[fv{1} fv{2} fv{3} fv{4}]=dwt2(cH3,'db9');
[fv{5} fv{6} fv{7} fv{8}]=dwt2(cV3,'db9');
[fv{9} fv{10} fv{11} fv{12}]=dwt2(cD3,'db9');
%%Level 2
[fv{13} fv{14} fv{15} fv{16}]=dwt2(cH2,'db9');
[fv{17} fv{18} fv{19} fv{20}]=dwt2(cV2,'db9');
[fv{21} fv{22} fv{23} fv{24}]=dwt2(cD2,'db9');
%%Level 1
[fv{25} fv{26} fv{27} fv{28}]=dwt2(cH1,'db9');
[fv{29} fv{30} fv{31} fv{32}]=dwt2(cV1,'db9');
[fv{33} fv{34} fv{35} fv{36}]=dwt2(cD1,'db9');
%%%Standard Deviation
for i=1:36
      fv{i}=std(fv{i}(:));
end
 if (exist('fp_database.dat')==2)
            load('fp_database.dat','-mat');
            fp_number=fp_number+1;
            data{fp_number,1}=fv;
            save('fp_database.dat','data','fp_number','-append');
        else
            fp_number=1;
            data{fp_number,1}=fv;
            save('fp_database.dat','data','fp_number');
 end
 msgbox('Image has been added to Database');
end

end


function hDelitemCallback(hObject, eventdata, handles)
    if (exist('fp_database.dat')==2)
            button = questdlg('Do you really want to remove the Database?');
            if strcmp(button,'Yes')
                delete('fp_database.dat');
                msgbox('Database was succesfully removed from the current directory.','Database removed','help');
            end
        else
            warndlg('Database is empty.',' Warning ')
    end
    delete('Database/*');
end

function hMatchitemCallback(hObject, eventdata, handles)
fingerprint=getappdata(hOpenMenuitem,'image');    
[XofCenter,YofCenter] = centralizing(fingerprint);
i=XofCenter-31;
j=YofCenter-31;
img=imcrop(fingerprint,[i j 63 63]);
fv=cell(1,36);
[cA1,cH1,cV1,cD1] = dwt2(img,'db9');
[cA2,cH2,cV2,cD2] = dwt2(cA1,'db9');
[cA3,cH3,cV3,cD3] = dwt2(cA2,'db9');
%%%extract features
%%Level 3
[fv{1} fv{2} fv{3} fv{4}]=dwt2(cH3,'db9');
[fv{5} fv{6} fv{7} fv{8}]=dwt2(cV3,'db9');
[fv{9} fv{10} fv{11} fv{12}]=dwt2(cD3,'db9');
%%Level 2
[fv{13} fv{14} fv{15} fv{16}]=dwt2(cH2,'db9');
[fv{17} fv{18} fv{19} fv{20}]=dwt2(cV2,'db9');
[fv{21} fv{22} fv{23} fv{24}]=dwt2(cD2,'db9');
%%Level 1
[fv{25} fv{26} fv{27} fv{28}]=dwt2(cH1,'db9');
[fv{29} fv{30} fv{31} fv{32}]=dwt2(cV1,'db9');
[fv{33} fv{34} fv{35} fv{36}]=dwt2(cD1,'db9');
%%%Standard Deviation
for i=1:36
      fv{i}=std(fv{i}(:));
end  
% Checking with DataBase
if (exist('fp_database.dat')==2)
    a=load('fp_database.dat','-mat');
    %whos fv          
    % start checking ---------------------------------------
    
    flag=0;
    for scanning=1:a.fp_number
    d=0;    
    fv1=a.data{scanning,1};
    fvt=cell2mat(fv);
    fv1t=cell2mat(fv1);
    %whos fv1t
    for i=1:36
        d=d+(fvt(i)-fv1t(i))*(fvt(i)-fv1t(i));
    end
    d=sqrt(d)
        if(d==0)
            msgbox('The Fingerprint is present in DataBase');
            flag=1;
            break;
        end
        
    end
    if(flag==0)	
        msgbox('The Fingerprint is not present in DataBase');
    end
else
            message='DataBase is empty. No check is possible.';
            msgbox(message,'FingerCode DataBase Error','warn');
end

end

function hAboutitemCallback(hObject, eventdata, handles)
msgbox('This project is done by Manit Kumar and Chidananda H.G from NMIT, Bangalore under the guidance of Dr. H. Sarojadevi, Professor of department of CSE, NMIT-Bangalore and Ms.Vidyadevi G Biradar, Assistant Professor of department of CSE, NMIT-Bangalore');
end


function hContactitemCallback(hObject, eventdata, handles)
msgbox('For any details contact Manit Kumar, Email:vikky_manit@yahoo.co.in or Chidananda H G, Email:mrchidu@gmail.com');
end

end


