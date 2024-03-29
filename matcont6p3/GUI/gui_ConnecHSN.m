function gui_ConnecHSN

global gds
    gds.gui_point = @point;
    gds.gui_type = @curvetype;
    gds.gui_starter = @starter1;
    gds.gui_load_layout = @load_layout;
    gds.gui_layoutbox = @layoutbox;
    gds.gui_numeric = @numeric1;
    gds.gui_load_draw = @load_draw;
    gds.gui_load_attributes = @load_attributes;
    gds.gui_make_ready = @make_ready;
    gds.gui_label = @label1;
    gds.gui_numeric_label = @numeric_label;
    gds.gui_draw = @draw;
    gds.gui_output = @output;
    gds.gui_load_point = @load_point;
    gds.gui_singularities = @singularities;
    gds.gui_start_cont = @start_cont;

%-------------------------------------------------------------------------
function point(tag)       
%dit roep je op als je zegt, initial point, point, en dan als Curve Con pakt!
global gds

[point,str]=strtok(tag,'_');
set(0,'ShowHiddenHandles','on');
h=findobj('Type','uimenu','Tag','window');
set(h,'enable','on');
h=findobj('Type','uimenu','Tag','curvetype');
set(0,'ShowHiddenHandles','off');
list = get(h,'children');
type = get(list,'tag');
str  = strcat(str,'_');
for i=1:length(list)
    if isempty(findstr(str,type{i}))    
        set(list(i),'enable','off');        
    else
        set(list(i),'enable','on');       
    end
end
gds.point=sprintf('%s%c',point,char(32));
curvetype;
%-------------------------------------------------------------------------    
function curvetype
global gds path_sys MC

if isempty(gds.point)
    gds.point = 'P ';   
end
gds.type = 'ConnecHSN ';
matcont('make_curve');
load_matcont;
file = strcat(gds.system,'.mat');
file = fullfile(path_sys,file);
save(file,'gds');
integrator; starter;
if  ~isempty(MC.numeric_fig),numeric;end
set(0,'ShowHiddenHandles','on');        
h = findobj('Type','uimenu','Tag','window');
set(h,'enable','on');
set(0,'ShowHiddenHandles','off');        
%-------------------------------------------------------------------------    
function starter1(handles)
global gds HTHSNds

ndim = size(gds.parameters,1);
color = [1 1 1];
HTHSNds.index = 0;
gds.options.ActiveParams=[];
gds.options.ActiveUParams=[];
gds.options.ActiveSParams=[];
s = gds.dim*2+ndim*2+40;
slider_step(1) = 2/s;
slider_step(2) = 2/s;
set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
j = s-2;
stat = uicontrol(handles.figuur,'Style','text','String','Initial Point','Tag','Initial_Point','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
pos = [10 j 35 1.80]; user.num = 0; user.pos = pos;
set(stat,'Position',pos,'UserData',user);
pos1 = starter('start_time',handles,j);
pos1 = start_coordinatesx0(handles,pos1);
pos1 = starter('start_parameters_orbit',handles,pos1); 
pos1 = start_c(handles,pos1);
pos1 = start_period(handles,pos1);
pos1 = start_select_cycle(handles,pos1);%zet die pushbutton van select cycle erbij in het starter window
starter('in_start',handles);
%-------------------------------------------------------------------------    
function load_layout(handles)
global gds 

for i=1:3
    if gds.numeric.ConnecHSN{i,2}==1    
        gds.numeric.ConnecHSN{i,1}=upper(gds.numeric.ConnecHSN{i,1});        
    elseif gds.numeric.ConnecHSN{i,2}==0    
        gds.numeric.ConnecHSN{i,1}=lower(gds.numeric.ConnecHSN{i,1});        
    end
    string{i,1}=gds.numeric.ConnecHSN{i,1};   
end
set(handles.layoutbox,'String',string);
%-------------------------------------------------------------------------    
function layoutbox(list,index_selected)
global gds 

c = gds.numeric.ConnecHSN{index_selected,2};    
gds.numeric.ConnecHSN{index_selected,2} = 1-c;    
%-------------------------------------------------------------------------    
function numeric1(handles)
global gds
ndim=size(gds.parameters,1);

if ~isfield(gds.numeric,'ConnecHSN') 
    gds.numeric.ConnecHSN = {'time' 1;'coordinates' 1;'eps1' 1;'parameters' 0'};    
end

s = 2;
if gds.numeric.ConnecHSN{1,2}==1
    s = s+3;   
end
if gds.numeric.ConnecHSN{2,2}==1
    s = s+2*gds.dim+2;   
end
if gds.numeric.ConnecHSN{3,2}==1
    s = s+40;   
end
if gds.numeric.ConnecHSN{4,2}==1
    s = s+2*ndim+2;    
end
s = s+2;
slider_step(1) = 2/s;
slider_step(2) = 2/s;
set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
j=s;  
if gds.numeric.ConnecHSN{1,2}==1
    j=numeric('start_time',handles,j);    
end
if gds.numeric.ConnecHSN{2,2}==1
    j=numeric('start_coordinates',handles,j);    
end
if gds.numeric.ConnecHSN{3,2}==1
    j = numeric_eps1(handles,j);   
end
if gds.numeric.ConnecHSN{4,2}==1
    j=numeric('start_parameters',handles,j);    
end

numeric('in_numeric',handles);
 %-------------------------------------------------------------------------    
function [string,plotopts] = load_draw
global gds

plotopts = {};
string = {'Coordinates';'Parameters';'Time'};
if ~isempty(gds.options.Userfunctions) && gds.options.Userfunctions ~= 0
    string{end+1} = 'Userfunction';    
end
%-------------------------------------------------------------------------    
    
function load_attributes

%-------------------------------------------------------------------------


function e = make_ready(plo,x)

len=size(plo,2);
e(1,len) = 0;
for j=1:len
    switch plo(j).type   
        case 'Coordinates'        
            e(1,j)=plo(j).val;
        case 'Parameters'
            e(1,j)=-plo(j).val;    
        case 'Time'
            e(1,j)=0;
        case 'eps1'
            e(1,j)=1000;                
        case 'Userfunction'
            e(1,j) = -100000*ndim - pl(j).val;
        otherwise
            for j = 1:len            
                e(1,j)=inf;                
                return                                
            end
    end
end    
%-------------------------------------------------------------------------
function label = numeric_label(numsing) 
global gds MC

num = 1;
label{gds.dim} = '';    
if isfield(MC.numeric_handles,'coord1')
    for d=1:gds.dim    
        label{num}=sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(x(%d,i),''%%.8g''))',d,d);        
        num = num+1;       
    end
end
if isfield(MC.numeric_handles,'time')
    label{num}=sprintf('set(MC.numeric_handles.time,''String'',num2str(t(i,1),''%%.10g''))');     
    num = num+1;    
end
if isfield(MC.numeric_handles,'eps1')
    label{num}=sprintf('set(MC.numeric_handles.eps1,''String'',num2str(eps1(i,1),''%%.10g''))');     
    num = num+1;    
end
if isfield(MC.numeric_handles,'param1')
    for d=1:size(gds.parameters,1);    
        label{num}=sprintf('set(MC.numeric_handles.param%d,''String'',gds.parameters{%d,2})',d,d);        
        num = num+1;       
    end
end
%-------------------------------------------------------------------------
function draw(varargin)
global gds MC

ndim = size(gds.parameters,1);

if strcmp(varargin{4},'select')      
    if nargin < 7    
        x = varargin{3}; s = varargin{5};t = varargin{2};         
    else
        x = varargin{3}; s = varargin{5};h = varargin{6}; f = varargin{7};        
    end
else
    file = varargin{3};load(file);    
    s(1)   = [];   s(end) = [];   
end
switch varargin{1}
    case 2    
        plot2=varargin{2};        
        if (plot2==[inf inf]),return;end        
        hold on;d=axis;        
        p=0;axis(d);
        plo2{1}='empt';plo2{2}='empt' ;
        k=size(x,2);
        plo2s{1} = 'empt';plo2s{2} = 'empt';      
        p = size(cat(1,s.index),1);
        if p>0,        sdata = cat(1,s.data); end

        for j=1:2
            
            if plot2(1,j) <= -100000*ndim % case Userfunction
                plo2{j}  = h(2-plot2(1,j)-100000*ndim,:);
            
            elseif plot2(1,j)<0
                pl2=-plot2(1,j);
                plo2{j}=gds.parameters{pl2,2}*ones(1,k);
                if p>0, plo2s{j} = gds.parameters{pl2,2}*ones(1,p);end
            end      
            if plot2(1,j)==0
                plo2{j}=t(1:k);
                if p>0, plo2s{j}=cat(1,sdata.t);end
            end
            if strcmp(plo2{j},'empt')==1
                plo2{j}=x(plot2(1,j),1:k);
                if p>0
                    sdatax=cat(1,sdata.x)';
                    plo2s{j} =  sdatax(plot2(1,j),:);
                end
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo2{1},'empt')&&~strcmp(plo2{2},'empt')
              line(plo2{1},plo2{2},'LineStyle','-','Color','b');
            end
           if ~strcmp(plo2s{1},'empt')&&~strcmp(plo2s{2},'empt') & p>0
               line(plo2s{1},plo2s{2},'LineStyle','none', 'Marker' ,'*','Color','r');
           end     
        else
            plot(plo2s{1},plo2s{2},'m+','MarkerSize',200);                
             hold off;
        end
    case 3 %3D
        plo=varargin{2};
        if (plo==[inf inf inf]),return;end
        hold on;d=axis;
        p=0;axis(d);
        skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
        k = size(x,2);p=size(cat(1,s.index),1);
         if p>0,        sdata = cat(1,s.data);end
        axis(d);plo3{1} = 'empt';plo3{2} = 'empt';plo3{3} = 'empt';
        plo3s{1} = 'empt'; plo3s{2} = 'empt';plo3s{3} = 'empt';
        k=size(x,2);
        for j=1:3
            if plo(1,j) <= -100000*ndim % case Userfunction
                plo3{j}  = h(2-plo(1,j)-100000*ndim,:);
            
            elseif plo(1,j)<0
                pl3=-plo(1,j);
                plo3{j}  = gds.parameters{pl3,2}*ones(1,k);
                if p>0, plo3s{j} = gds.parameters{pl2,2}*ones(1,p);end
            end               
            if plo(1,j)==0
                plo3{j}  = t(1:k);
                if p>0,plo3s{j} = cat(1,sdata.t);end
            end
            if strcmp(plo3{j},'empt')==1
                plo3{j}  = x(plo(1,j),1:k);
                if p>0
                    sdatax   = cat(1,sdata.x)';
                    plo3s{j} =  sdatax(plo(1,j),:);
                end
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo3{1},'empt')&&~strcmp(plo3{2},'empt')&&~strcmp(plo3{3},'empt')
                line(plo3{1},plo3{2},plo3{3},'linestyle','-','Color','b');
            end
             if ~strcmp(plo3s{1},'empt')&&~strcmp(plo3s{2},'empt')&&~strcmp(plo3s{3},'empt')& p>0
                line(plo3s{1},plo3s{2},plo3s{3},'linestyle','none','Marker','*','Color','r');
             end
        else
            plot3(plo3s{1},plo3s{2},plo3s{3},'m+','MarkerSize',200);                
            hold off;
        end
    case 4 %numeric
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        dat = get(MC.numeric_fig,'Userdata');
        i = s.index;
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end
end

%-------------------------------------------------------------------------
function output(numsing,xout,s, hout,fout,i)
global MC gds
if ~isempty(MC.numeric_fig) %numeric window is open
    ndim = size(gds.parameters,1);    
    parameters = zeros(1,ndim);    
    for d = 1:ndim    
        parameters(d) = gds.parameters{d,2};           
    end                        
    
    dat = get(MC.numeric_fig,'Userdata');
    for k=1:size(dat.label,2)    
        eval(dat.label{k});      
    end
end
%-------------------------------------------------------------------------    
function load_point(index,x,string,file,varargin)
global gds HTHSNds

load(file);

for i = 1:gds.dim
    gds.coordinates{i,2} = x(i,index);   
end
gds.discretization.ntst = HTHSNds.ntst;
gds.discretization.ncol = HTHSNds.ncol;
ndim = size(gds.parameters,1);
for i = 1:ndim
    gds.parameters{i,2} = HTHSNds.P0(i,1);    
end

gds.T = HTHSNds.T;
gds.eps0 = HTHSNds.eps0;
gds.eps1 = HTHSNds.eps1;
gds.x0 = HTHSNds.x0;

gds.UParams = HTHSNds.UParams;
gds.SParams = HTHSNds.SParams;
gds.extravec = HTHSNds.extravec;
      
if strcmp(string,'save')
    num = varargin{1};  
    if ~exist('v')
        v = [];        
    end
    if ~exist('s')
        s = [];
    end
    if ~exist('h')
        h = [];
    end
    if ~exist('f')    
        f = [];        
    end
    if ~exist('cds')    
        cds = [];        
    end
    if ~exist('ctype')    
        ctype = [];        
    end
    if ~exist('point')    
        point = [];        
    end
    save(file,'x','v','s','h','f','num','cds','HTHSNds','ctype','point');    
end    

function singularities  

%-------------------------------------------------------------------------
function start_cont(varargin)
global gds MC path_sys t dim_npos QS1 eig0 %% je gaat connBds ier ook toevoegen en hier de waarde toekennen aan de c's en tau's

func_handles = feval(gds.system);

if gds.c{1,2}==0
    error('UParam1 has to be different from zero');
end

if gds.eps0 == 0
    error('eps0 has to be different from zero');
end

if gds.eps0 < 0
    error('eps0 has to be positive');
end

check = 0;
if ~isempty(func_handles{3})
    jaco = func_handles{3};
    check = 1;
end
par = vertcat(gds.parameters{:,2});
p = num2cell(par);

if check ==1
    A = feval(jaco, 0, vertcat(gds.x0{:,2}), p{:});
else
    point = vertcat(gds.x0{:,2});
    x1 = point;
    x2 = point;
    for i=1:gds.dim
        x1(i) = x1(i)-gds.options.Increment;
        x2(i) = x2(i)+gds.options.Increment;
        j(:,i) = feval(func_handles{2}, 0, x2, p{:})-feval(func_handles{2}, 0, x1, p{:});
        x1(i) = point(i);
        x2(i) = point(i);
    end
    A = j/(2*gds.options.Increment);
end

[V,D] = eig(A); %dimensie van D : dim maal 1
[K,i] = min(abs(diag(D)));%eigenwaarde 0
eig0 = V(:,i); 

% is x_0 an equilibrium?
func = feval(func_handles{2}, 0, vertcat(gds.x0{:,2}), p{:});

if ~(norm(func)<=1e-1)
    errordlg('The starting point is not an equilibrium','Bad startpoint');   
    return
end

if ~(abs(K)<=1e-1)
    errordlg('The starting point is not a saddle-node','Bad startpoint');   
    return
end    

B = D([1:i-1 i+1:end],[1:i-1 i+1:end]);%je sluit eigenwaarde 0 uit
dim_nneg = sum(real(diag(B)) < 0);
if (dim_nneg == gds.dim-1)   
    if min(abs(real(diag(B)))) < 1e-10
        dim_nneg = dim_nneg -1;    
    end
end
if (dim_nneg == 0)
    if min(abs(real(diag(B)))) < 1e-10
        dim_nneg = dim_nneg +1;
    end
end
dim_npos = gds.dim-dim_nneg-1; 

gds.c{1,2} = gds.c{1,2}/abs(gds.c{1,2});
start_point = vertcat(gds.x0{:,2}) + gds.eps0*gds.c{1,2}*eig0;

ndim=size(gds.parameters,1);
set(0,'ShowHiddenHandles','on');
c=findobj('Style','edit','Tag','interval');
d=findobj('Style','edit','Tag','edit0');
set(0,'ShowHiddenHandles','off');
val=str2double(get(c,'String'));        
gds.time{1,2} = str2double(get(d,'String'));    
gds.integrator.tspan(2) = gds.time{1,2}+val;
gds.integrator.tspan(1) = gds.time{1,2};

if ~isempty(varargin)     
    funhandle = feval(gds.system);    
    file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));    
    if exist(file,'file')    
        load(file);    
    else
        errordlg('It is not possible to extend the current curve!','Error extend P');        
        return 
        
    end
    x0 = x(1:gds.dim,end);    
    for k = 1:ndim   
        gds.parameters{k,2} = param(k);       
    end
    if exist('option','var')    
        gds.integrator.options = option;
        
    end
    set(0,'ShowHiddenHandles','on');    
    hh = findobj('Style','edit','tag','interval');   
    val = str2double(get(hh,'String'));    
    gds.time{1,2} = abs(t(end));    
    if val<0 then    
        errordlg('this isn''t possible');        
        set(hh,'string','1');        
    else
        gds.integrator.tspan=[gds.time{1,2} (val+gds.time{1,2})];        
    end
else
    funhandle = feval(gds.system);    
    x0 = start_point; %hij moet van dat start punt beginnen   
    gds.integrator.options = odeset(gds.integrator.options,'Jacobian',funhandle{3});    
    gds.integrator.options = odeset(gds.integrator.options,'Hessians',funhandle{5});            
end

try
    if (gds.options.Backward==1)   
        interval = abs(gds.integrator.tspan(2)-gds.integrator.tspan(1));        
        tspan = [gds.integrator.tspan(1) gds.integrator.tspan(1)-interval];        
    else
        tspan = gds.integrator.tspan;        
    end
    str = sprintf('(funhandle{2},tspan,x0,gds.integrator.options');   
    par = cellstr(strcat(',',num2str(vertcat(gds.parameters{:,2}))));    
    str = strcat(str,strcat(par{:,1}),')');      
    
    funhandle = feval(gds.system);
    
    switch gds.integrator.method    
        case 'ode45'                      
            str = strcat('ode45',str);
        case 'ode23'
            str = strcat('ode23',str);
        case 'ode113'
            str = strcat('ode113',str);
        case 'ode15s'
            str = strcat('ode15s',str);
        case 'ode23s'
            str = strcat('ode23s',str);
        case 'ode23t'
            str = strcat('ode23t',str);
        case 'ode23tb'
            str = strcat('ode23tb',str);
        case 'ode78'
            str = strcat('ode78',str);
        case 'ode87'
            str = strcat('ode87',str);        
    end
    npoints = 1;    
    %%%    
    if size(MC.D2,2)>0||size(MC.D3,2)>0||~isempty(MC.numeric_fig)    
        gds.integrator.options = odeset(gds.integrator.options,'OutputFcn',@integ_plotConnecHSN);        
    else
        gds.integrator.options = odeset(gds.integrator.options,'OutputFcn',@integ_prs);        
    end

    file = fullfile(path_sys,gds.system);    
    save(file,'gds');    
    set(0,'ShowHiddenHandles','on');    
    duration = findobj('style','text','Tag','duration');    
    status = findobj('style','text','tag','mstatus');    
    set(0,'ShowHiddenHandles','off');    
    set(status,'string','computing');    
    StartTime = clock; %global t
           
    [t,y] = eval(str);      
            
    y1 = zeros(size(y,1),gds.dim+1);   
    for k = 1:size(y,1)    
        y1(k,1:gds.dim) = y(k,:);        
        y1(k,gds.dim+1) = norm(vertcat(gds.x0{:,2})'-y(k,:));        
    end

    EndTime = clock;    
    string  = sprintf('%.1f secs\n', etime(EndTime, StartTime));    
    set(duration,'String',string);    
    set(status,'String','ready');    
    x = y1';    
    file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);    
    if isempty(varargin)    
        s(1,1).index = 1;         s(1,1).msg = 'This is the first point of the orbit'; s(1,1).label = '00';       
        s(1,1).data.t = t(1); s(1,1).data.x = x(1:gds.dim,1)';        
        s(2,1).index = size(t,1); s(2,1).msg = 'This is the last point of the orbit';  s(2,1).label = '99';        
        s(2,1).data.t = t(end); s(2,1).data.x = x(1:gds.dim,end)';        
    else %dus als je extend doet    
        x1 = x; t1 = t;        
        load(file);        
        x = [x,x1];        
        t = [t;t1];        
        s(1,1).index = 1;         s(1,1).msg = 'This is the first point of the orbit'; s(1,1).label = '00';        
        s(2,1).index = size(t,1); s(2,1).msg = 'This is the last point of the orbit';  s(2,1).label = '99';        
        s(2,1).data.t = t(end); s(2,1).data.x = x(1:gds.dim,end)';       
    end

    [QS1, eigvl1] = computeBaseConnecHSN(A);%n maal n matrix
    
    ctype  = gds.type; point = gds.point;            
    status = mkdir(path_sys,gds.system);    
    dir = fullfile(path_sys,gds.system);    
    status = mkdir(dir,gds.diagram);    
    param  = vertcat(gds.parameters{:,2});      
    option = gds.integrator.options;    
    save(file,'t','x','s','param','ctype','point','option');     
    file=fullfile(path_sys,gds.system);save(file,'gds');    
    sys.gui.plot_points = npoints;   
catch
    errordlg(lasterr,'error integration');    
end 
    
%-----------------------------------------------------------------------------------    
function j = numeric_eps1(handles,j)
global gds
color = [1 1 1];
j = j-2;user.num=0;
stat = uicontrol(handles.numeric_fig,'Style','text','String','eps1','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [12 j 18 1.8];user.pos=pos;
set(stat,'Position',pos,'UserData',user);

stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','eps1','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'units','characters');
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','eps1','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'units','characters');
j = j-2;
pos = [2 j 18 1.8];user.pos=pos;
set(stat1,'Position',pos,'UserData',user);
pos = [20 j 27 1.8];user.pos=pos;
set(stat2,'Position',pos,'UserData',user);
guidata(handles.numeric_fig,handles);

%-----------------------------------------------------------------------------------
function  [Q0,sortedevls] =  computeBaseConnecHSN(A0,flag)
%unstable
        [evc0,evl0] = eig(A0);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        
        evl = diag(evl);
        pos = find(real(evl)<0);
        evlsr = real(evl(pos));
        evls = evl(pos);
        evcs = evc(:,pos); % the eigenvectors involved
        
        [a,b] = sort(evlsr); % result: a(i) = val(b(i))
        VU = evcs(:,b);
        sortedevls = evls(b);                    
    
        VU_1 = zeros(size(A0,1),size(VU,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(VU,2)
            VU_1(:,f) = VU(:,end-f+1);
            sortedevls_1(f) = sortedevls(end-f+1);
        end
        sortedevls = sortedevls_1';
    
        B = VU_1;
        k = 1;
        while k <= size(B,2)
            ok = 1;
            init = 1;
            while ok == 1 && init <= size(B,1)
                if imag(B(init,k))                     
                    tmp = B(:,k);
                    B(:,k) = real(tmp);
                    B(:,k+1) = imag(tmp);
                    ok = 0;
                end
                init = init+1;
            end
            if ok == 1
               k = k+1; 
            else
                k = k+2;
            end            
        end
        
        % Compute orthonormal basis for the eigenspace
        [Q0,RU] = qr(B);
    
%--------------------------------------------------------------------------
function j = start_c(handles,j)
global gds path_sys;
color = [1 1 1];
user.num = 0;
gds.c = [];
gds.c{1,1} = 'UParam1';
gds.c{1,2} = 0;


   if strcmp(gds.c{1,2},'')
       string = '';
   else
       string = num2str(gds.c{1,2});
   end
   tag1 = strcat('text',num2str(1));
   tag2 = strcat('cc',num2str(1));
   
   j = j-2;
   pos = [2 j 18 1.8];
   user.pos = pos;
   stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.c{1,1},'Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12,'UserData',user);
   set(stat,'Position',pos,'UserData',user);
      
   pos = [20 j 25 1.8];user.pos=pos;
   edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
   set(edit,'Callback','Callback');         
   set(edit,'Position',pos,'UserData',user);

guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j = start_coordinatesx0(handles,j)
%enters the names of the coordinates in the window
%it the value is known, it is entered
global gds;
color = [1 1 1];
user.num = 0;
ndim = size(gds.parameters,1);
gds.x0 = [];
for i = 1:gds.dim
    gds.x0{i,1} = strcat('x0_',gds.coordinates{i,1});
    gds.x0{i,2} = 0;
end

for i = 1:(gds.dim)
    if strcmp(gds.coordinates{i,1},'')
        string = '';
    else
        string = num2str(gds.coordinates{i,2},'%0.8g');
        gds.x0{i,2} = gds.coordinates{i,2};
    end
    tag1 = strcat('text',num2str(gds.dim+ndim+i));
    tag2 = strcat('edit',num2str(gds.dim+ndim+i));
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','Callback');
    j = j-2;
    pos = [2 j 18 1.8];
    user.pos = pos;
    stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.x0{i,1},'Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(stat,'Position',pos,'UserData',user);
    pos = [20 j 25 1.8];user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
end

guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j = start_select_cycle(handles,j)
%enter field for select cycle button
global gds path_sys;
j = j-3;
pos  = [10 j 30 2]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','pushbutton','String','Select Connection','Tag','select_cycleConnecHSN','units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user,'Callback','selectcycleConnecHSN');
guidata(handles.figuur,handles);
%--------------------------------------------------------------------------
function j = start_period(handles,j)
global gds;
gds.eps0 = 0;

color = [1 1 1];
j = j-2;
post = [2 j 18 1.8]; user.num = 0; user.pos = post; pose = [20 j 25 1.8];
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','eps0','Tag','eps0','Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12,'UserData',user);
set(stat,'Position',post,'UserData',user);
user.pos = pose;
edit  = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',num2str(gds.eps0),'Tag','eps0val','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','Callback');
set(edit,'Position',pose,'UserData',user);

guidata(handles.figuur,handles);
