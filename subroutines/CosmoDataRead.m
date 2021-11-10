function [num,sampName,X,DEMdata] = CosmoDataRead(filename,varargin)

[~,sheets] = xlsfinfo(filename);

if ~isempty(varargin)
    ind  = varargin{1}; 
else
    ind  = 1;
end

nuclide_sheet = ~cellfun(@isempty,strfind(sheets, 'Cronus'));    % find the excel sheet that have nuclide information

switch length(sheets)
    case 2           % single CRN plus Weathering
        nuclide = sheets{nuclide_sheet}(1:4);                            % identify which nuclide was provided
        switch nuclide, case '10Be', X.n = 1; case '36Cl', X.n = 2; end  % save which nuclide is being used  
        
        % load data
        [num,txt,~]    = xlsread(filename,find(nuclide_sheet));          % load CRN data
        [Xdata,Xtxt,~] = xlsread(filename,find(~nuclide_sheet));         % load compositional data,
        
        X.W              = Xdata(10);                                    % Weathering rate
        X.Wstd           = Xdata(11);
        
        num = num(ind,:); 
        sampName = txt{ind+2,1}; 
        
    case 3          % paired CRN
        [num10,txt10,~] = xlsread(filename,'10Be Cronus');               % load 10Be data
        [num36,~,~]     = xlsread(filename,'36Cl Cronus');               % load 36Cl data
        [Xdata,Xtxt,~]  = xlsread(filename,find(~nuclide_sheet));        % load compositional data
        
        num10 = num10(ind,:);   
        num36 = num36(ind,:);
        sampName = txt10{ind+2,:}; 
        num.num10 = num10; num.num36 = num36;
end

Xdata = Xdata(ind,:);                                                    % select sample data
DEMdata.method = Xtxt{ind+1,10};
X.soil_mass    = Xdata(7);
global scaling_model
scaling_model  = Xtxt{ind+1,9};       


if isnan(Xdata(2))
    X.fQzS = Xdata(1); 
    X.fCaS = Xdata(2); 
    X.fXS  = Xdata(3);  
    X.mode = 'soil';
else
    X.fQzB = Xdata(1); 
    X.fCaB = Xdata(2); 
    X.fXB  = Xdata(3);  
    X.mode = 'bedrock';
end


end

