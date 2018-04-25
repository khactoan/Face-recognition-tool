function varargout = FaceRecognitionTool(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FaceRecognitionTool_OpeningFcn, ...
    'gui_OutputFcn',  @FaceRecognitionTool_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function FaceRecognitionTool_OpeningFcn(hObject, eventdata, handles, varargin)
axes(handles.axes4)
cla
axes(handles.axes3)
cla
handles.output = hObject;
guidata(hObject, handles);

function varargout = FaceRecognitionTool_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function pushbutton4_Callback(hObject, eventdata, handles)
global A m1 n1 No_Files_In_Class_Folder Class_Count Training_Set_Folder

Training_Set_Folder = [uigetdir(''),'\'];
m1 = 8;
n1 = 4;
TS_Vector = dir(Training_Set_Folder);
No_Folders_In_Training_Set_Folder = length(TS_Vector);
File_Count = 1;
Class_Count = 1;
h = waitbar(0,'Reading Test Images,Please wait...');
for k = 3:No_Folders_In_Training_Set_Folder
    waitbar(k/(No_Folders_In_Training_Set_Folder-2))
    Class_Folder = [Training_Set_Folder '\' TS_Vector(k).name,'\'];
    CF_Tensor = dir(Class_Folder);
    No_Files_In_Class_Folder(Class_Count) = length(CF_Tensor)-2;
    %     strr = sprintf('Reading Test Images...!, # of Classes = %d, Now Reading %d ',No_Folders_In_Training_Set_Folder-2,Class_Count);
    %     set(handles.edit3,'String',strr);
    drawnow;
    for p = 3:No_Files_In_Class_Folder(Class_Count)+2
        Tmp_Image_Path = Class_Folder;
        Tmp_Image_Name = CF_Tensor(p).name;
        Tmp_Image_Path_Name = [Tmp_Image_Path,Tmp_Image_Name];
        if strcmp(Tmp_Image_Name,'Thumbs.db')
            break
        end
        test = imread(Tmp_Image_Path_Name);
        if length(size(test))==3
            Tmp_Image = rgb2gray(test);
        else
            Tmp_Image = test;
        end
        Tmp_Image_Down_Sampled = double(imresize(Tmp_Image,[m1 n1]));
        Image_Data_Matrix(:,File_Count) = Tmp_Image_Down_Sampled(:);
        File_Count = File_Count+1;
        
    end
    Class_Count = Class_Count+1;
    
end
close(h)
A = Image_Data_Matrix;
A = A/(diag(sqrt(diag(A'*A))));


function pushbutton5_Callback(hObject, eventdata, handles)
global A m1 n1 No_Files_In_Class_Folder Class_Count Training_Set_Folder
d=uigetdir('','Select Input-folder'); %select the input-folder that contains the subfolders
cd(d);
list = dir;
list = list([list.isdir]); 
list = list(~ismember({list.name},{'.' '..'}));
listSize=length(list);
global countDetected;
countDetected = 0;
tic;
for folder=1:listSize
      oldfolder = cd(list(folder).name);
      files = dir('*.pgm');
      numberOfFiles = length(files);
      for k = 1:numberOfFiles
          file_path = strcat(pwd, '\', files(k).name);
          test = imread(file_path);
          if length(size(test))==3
                Test_Image = rgb2gray(test);
            else
                Test_Image = test;
            end



            Test_Image_Down_Sampled = double(imresize(Test_Image,[m1 n1]));
            y = Test_Image_Down_Sampled(:);
            n = size(A,2);
            f=ones(2*n,1);
            Aeq=[A -A];
            lb=zeros(2*n,1);
            x1 = linprog(f,[],[],Aeq,y,lb,[],[],[]);
            x1 = x1(1:n)-x1(n+1:2*n);
            nn = No_Files_In_Class_Folder;
            nn = cumsum(nn);
            tmp_var = 0;
            k1 = Class_Count-1;
            for i = 1:k1
                delta_xi = zeros(length(x1),1);
                if i == 1
                    delta_xi(1:nn(i)) = x1(1:nn(i));
                else
                    tmp_var = tmp_var + nn(i-1);
                    begs = nn(i-1)+1;
                    ends = nn(i);
                    delta_xi(begs:ends) = x1(begs:ends);
                end
                tmp(i) = norm(y-A*delta_xi,2);
                tmp1(i) = norm(delta_xi,1)/norm(x1,1);
            end
            Sparse_Conc_Index = (k1*max(tmp1)-1)/(k1-1);
            clss = find(tmp==min(tmp));
            cccc = dir([Training_Set_Folder]);
            trainFolder = cccc(clss+2).name;
            if trainFolder == list(folder).name
               countDetected = countDetected + 1; 
            else
               disp(file_path);
            end
      end
      cd(oldfolder); 
end
toc;
disp(['Recognition rate: ', num2str(countDetected/200 * 100), '%']);

function pushbutton6_Callback(hObject, eventdata, handles)
global A m1 n1 No_Files_In_Class_Folder Class_Count Training_Set_Folder
[Test_File Test_File_Path] = uigetfile('*.jpg;*.pgm;*.png;*.tif','Select a Test Image');
parentFolder = Test_File_Path(end-3:end-1);

test_image_path = [Test_File_Path Test_File];
axes(handles.axes3)
cla
axes(handles.axes4)
cla
axes(handles.axes3)
imshow(test_image_path);
title('Test Image','Color','red','FontSize',15);
drawnow;
Test_File = [Test_File_Path Test_File];
test = imread(Test_File);
if length(size(test))==3
    Test_Image = rgb2gray(test);
else
    Test_Image = test;
end
Test_Image_Down_Sampled = double(imresize(Test_Image,[m1 n1]));
y = Test_Image_Down_Sampled(:);
n = size(A,2);
f=ones(2*n,1);
Aeq=[A -A];
lb=zeros(2*n,1);
x1 = linprog(f,[],[],Aeq,y,lb,[],[],[]);
x1 = x1(1:n)-x1(n+1:2*n);
nn = No_Files_In_Class_Folder;
nn = cumsum(nn);
tmp_var = 0;
k1 = Class_Count-1;
for i = 1:k1
    delta_xi = zeros(length(x1),1);
    if i == 1
        delta_xi(1:nn(i)) = x1(1:nn(i));
    else
        tmp_var = tmp_var + nn(i-1);
        begs = nn(i-1)+1;
        ends = nn(i);
        delta_xi(begs:ends) = x1(begs:ends);
    end
    tmp(i) = norm(y-A*delta_xi,2);
    tmp1(i) = norm(delta_xi,1)/norm(x1,1);
end
Sparse_Conc_Index = (k1*max(tmp1)-1)/(k1-1);
clss = find(tmp==min(tmp));
cccc = dir([Training_Set_Folder]);
Which_Folder = dir([Training_Set_Folder,cccc(clss+2).name,'\']);
Which_Image = randsample(3:length(Which_Folder),1);
Image_Path = [Training_Set_Folder,cccc(clss+2).name,'\',Which_Folder(Which_Image).name];
Class_Image = (Image_Path);
axes(handles.axes4);
imshow(Class_Image)
if  cccc(clss+2).name == parentFolder
    title('Detected Succesfully','Color','blue','FontSize',15)
else
    title('Detect Failed','Color','red','FontSize',15)
end
