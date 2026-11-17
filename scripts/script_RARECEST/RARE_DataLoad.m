%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Analysis tool for CEST-Maps acquired with the cestRARE sequence
%%%
%%% Project:
%%%             Chemical Exchange Saturation Transfer
%%% 
%%% Author:
%%%             AB (V 1.0)
%%%             PB (V 1.1 and further)
%%% 
%%% Location:
%%%             Westfälische Wilhelms-Universität Münster
%%%             University Hospital Münster
%%%             Research Group for Experimental Magnetic Resonance
%%%
%%% Date:
%%%             V 1.0    February, 26th 2015
%%%                      --> project shifted from AB to PB
%%%             V 1.1    October, 20th 2015 [adaption of subprogramms for
%%%                      the calculation of corrected spectra; PB]
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function varargout = RARE_DataLoad(varargin)
% RARE_DATALOAD M-file for RARE_DataLoad.fig
%      RARE_DATALOAD, by itself, creates a new RARE_DATALOAD or raises the existing
%      singleton*.
%
%      H = RARE_DATALOAD returns the handle to a new RARE_DATALOAD or the handle to
%      the existing singleton*.
%
%      RARE_DATALOAD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RARE_DATALOAD.M with the given input arguments.
%
%      RARE_DATALOAD('Property','Value',...) creates a new RARE_DATALOAD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RARE_DataLoad_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RARE_DataLoad_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RARE_DataLoad

% Last Modified by GUIDE v2.5 26-Feb-2015 11:39:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RARE_DataLoad_OpeningFcn, ...
                   'gui_OutputFcn',  @RARE_DataLoad_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before RARE_DataLoad is made visible. 
function RARE_DataLoad_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RARE_DataLoad (see VARARGIN)

% Choose default command line output for RARE_DataLoad
handles.output = hObject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes 2016 01 21 - PB
% Creation of an automated file loading
%%%%%
% input('Bitte CEST Daten mit allen Offsets wählen [Enter drücken]')
% handles.CEST.CEST = read_paravision_data;
% input('Bitte CEST Daten ohne Offset [Enter drücken]')
% handles.CEST.CEST_S0 = read_paravision_data;
fprintf('%s \n','Please specify the scan IDs of the CEST measurements');
handles.pathname = uigetdir('*.*', 'Select scan ID folder');
fprintf('%s \n','Please specify the scan IDs of the CEST measurements');
handles.numbers.CEST_ao = input('Scan ID of the CEST measurement with all offsets: ');
handles.numbers.CEST_so = input('Scan ID of the CEST measurement with S_0 single offset acquisition: ');
handles.path_CEST_so = [handles.pathname '\' num2str(handles.numbers.CEST_so) '\pdata\1\2dseq'];
handles.path_CEST_ao = [handles.pathname '\' num2str(handles.numbers.CEST_ao) '\pdata\1\2dseq'];
handles.CEST.CEST = read_paravision_data(handles.path_CEST_ao);
handles.CEST.CEST_S0 = read_paravision_data(handles.path_CEST_so);
% End of changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
imagesc(squeeze(handles.CEST.CEST_S0.dataset(1,:,:)));
colormap('gray');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RARE_DataLoad wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RARE_DataLoad_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_ROI.
function pushbutton_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%% Steuern des ROI Buttons und erstellen der ROI Matrix %%%
handles.ROI.ROI_X = 0;
handles.ROI.ROI_Y = 0;

[handles.ROI.BI,...
    handles.ROI.ROI_X,...
    handles.ROI.ROI_Y] =roipoly;

handles.ROI.nonzero =...
    find(handles.ROI.BI);
handles.ROI.pixels = ...
    length(handles.ROI.nonzero);

handles.ROI.line = ...
    line(handles.ROI.ROI_X,...
    handles.ROI.ROI_Y, 'Color', 'red', 'LineWidth', 1);

handles.ROI.Data = ...
    zeros(handles.CEST.CEST.reco.RECO_size(1),handles.CEST.CEST.reco.RECO_size(2));

display('test');

%%% Update handles structure %%%
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% CEST Spektra
RARE_CEST_Spectra;
guidata(hObject, handles);



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% CEST Contrast
RARE_CEST_Contrast;
guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% CEST Contrast Corrected
RARE_CEST_Contrast_Corrected;
guidata(hObject, handles);


% --- Executes on button press in pushbutton5. Load WASSR Data
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RARE_CEST_LoadWASSRdata;
guidata(hObject, handles);


% --- Executes on button press in pushbutton6. CEST Contrast ROI
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RARE_CEST_Contrast_ROI;
guidata(hObject, handles);


% --- Executes on button press in pushbutton7. Contrast ROI Corrected
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RARE_CEST_Contrast_ROI_Corrected;
guidata(hObject, handles);


% --- Executes on button press in pushbutton8_spectra_corrected.
function pushbutton8_spectra_corrected_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_spectra_corrected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%SpectraCorrectedPixel;
RARE_SpectraCorrectedPixel;
guidata(hObject, handles);


% --- Executes on button press in pushbutton9_spectra_corrected_ROI.
function pushbutton9_spectra_corrected_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9_spectra_corrected_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%SpectraCorrectedROI;
RARE_SpectraCorrectedROI;
guidata(hObject, handles);

