function [correctedDistance,originalDistance,ellipseParam,signals] = ellipse_fit_single(SPD,CPD,MCPD,fit)
    % Function to do some 'simple' ellipse fitting on both the output lissajous
    % of HoQI. 
    % Syntax is [correctedDistance,originalDistance,ellipseParam,signals] =ellipse_fit_single(SPD,CPD,MCPD)
    % The script will sum the photodiode inputs and create a intensity
    % fluctuation suppressed lissajous. 

    if ~exist('MCPD','var') || isempty(MCPD)
        % Is the interferometer only running a single diode? if so we shouldn't
        % subtract something that doesn't exist, so we'll modify how we do the
        % atan2, lets set the MCPD to something obvious 
        MCPD=0;
    end

    if ~exist('fit','var') || isempty(fit)
        fit=struct();
        statusBad=1;                % Status is already bad so there's no point running validation
    else
        statusBad = checkFittingStruct(fit);
    end
    
    % Hard setting the wavelength and ensuring that the variables are the
    % same dimensions
    
    Lambda=1064e-9;
    [ms,ns] = size(SPD);
    if (ms < ns) 
       SPD=SPD.';
       CPD=CPD.';
       MCPD=MCPD.';
    end
    
    % Making a copy of the original signals
    SPD_orig = SPD;
    CPD_orig = CPD; 
    MCPD_orig =MCPD;
    

    if MCPD==1 
        % Here we're using three diodes so we'll use all three
        SPD_new = SPD-CPD;
        CPD_new = SPD-MCPD;
    else 
        % Only using two diodes so lets not confuse the fitting script.
        SPD_new = SPD;
        CPD_new = CPD;
    end
    
    if statusBad==1
        % Should not run the fitting script on the data below as the status of the supplied struct is bad.
        
        % Carrying out the original fitting 
       
        % Here we start the least squares fit
    
        X=[SPD_new.^2, SPD_new.*CPD_new, CPD_new.^2, SPD_new,CPD_new];      % Making a vector of the variables to fit

        Coeffs1= -sum(X)/(X'*X);                                            % calculating the coefficients

        [A1,B1,C1,D1,E1]=deal(Coeffs1(1),Coeffs1(2),Coeffs1(3),Coeffs1(4),Coeffs1(5));

        Alpha1=0.5*atan2(B1,(A1-C1));                                       % calculating the rotation
        a1=sqrt(abs(((A1+C1)-sqrt(B1.^2+(A1-C1).^2)))/2);                   % calculating the semi major axis
        b1=sqrt(abs(((A1+C1)+sqrt(B1.^2+(A1-C1).^2)))/2);                   % calculating the semi minor axis
        x01=(2*C1.*D1-B1.*E1)/(B1.^2-4*A1.*C1);                             % calculating the X offset
        y01=(2*A1.*E1-D1.*B1)/(B1.^2-4*A1.*C1);                             % calculating the Y offset


        K1=sqrt(sqrt(abs(A1.*(x01)^2+B1.*x01.*y01+C1.*(y01).^2-1))/((a1.*b1).^2));  % Calculating the scaling constants
        a1scale=K1*a1;                                                      % rescaling the parameters so they equal unity
        b1scale=K1*b1;                                                      % rescaling the parameters so they equal unity
       
    else
        
        % We can use the previous fit for our ellipse parameters
        a1scale=fit.a;
        b1scale=fit.b;
        Alpha1=fit.Rotation;
        x01=fit.xOffset;
        y01=fit.yOffset;
        K1=fit.ScalingConstant;
    end

    % Transforming the original parameters to reduce the non linearities
    SPD_cor=((SPD_new-x01).*cos(Alpha1)+(CPD_new-y01).*sin(Alpha1))./a1scale;
    CPD_cor=(-(SPD_new-x01).*sin(Alpha1)+(CPD_new-y01).*cos(Alpha1))./b1scale;
    
    SPD_cor = SPD_cor./max(SPD_cor);
    CPD_cor = CPD_cor./max(CPD_cor);
    ellipseParam = struct('a',a1scale,'b',b1scale,'Rotation',Alpha1,'xOffset',x01,'yOffset',y01,'ScalingConstant',K1);

    % Original signal calibration

    % Offsets
    x_off = (max(SPD_new)+min(SPD_new))/2;
    y_off = (max(CPD_new)+min(CPD_new))/2;

    % Correcting the gains
    x_gain = max(SPD_new)-min(SPD_new);
    y_gain = max(CPD_new)-min(CPD_new);

    % Now rescaling the original calibration signals
    SPD_cal = (SPD_new-x_off)./x_gain;
    CPD_cal = (CPD_new-y_off)./y_gain;
    
%     figure
%     plot(SPD_cal,CPD_cal)
    
    % Outputting the signals
    signals = struct();
    signals.fitted.SPD = SPD_cor;
    signals.fitted.CPD = CPD_cor;
    
    % Outputting the corrected (ellipes fitted distance) and the original
    % distance
    correctedDistance=(Lambda/(4*pi))*unwrap(atan2(SPD_cor,CPD_cor));
    originalDistance= (Lambda/(4*pi))*unwrap(atan2(SPD_cal,CPD_cal));

    function statusBad = checkFittingStruct(fit)
        % Function will check the validity of the fitting struct and ensure
        % all the paramaters are here. 
        statusBad=0; % We assume its ok. 
        if isfield(fit,'a')==0
            statusBad=1;
        end

        if isfield(fit,'b')==0
            statusBad=1;
        end

        if isfield(fit,'Rotation')==0
            statusBad=1;
        end

        if isfield(fit,'xOffset')==0
            statusBad=1;
        end

        if isfield(fit,'yOffset')==0
            statusBad=1;
        end

        if isfield(fit,'ScalingConstant')==0
            statusBad=1;
        end
        
    end

end