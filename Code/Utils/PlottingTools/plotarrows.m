function h = plotarrows(x,y,dx,dy,mx,my,scale,varargin)
%
% h = plotarrows(x,y,dx,dy,mx,my,scale,varargin)
%
% Plots data vectors, as well as, optionally, model vectors and data
% error ellipses.
% This code is intended to replace plot_map_vectors.m, arrows2.m,
% arrows_kyle.m, etc.
%
% By K. Anderson, spring 2012
%
% INPUT
%           x, y:       Station positions
%           dx, dy:     Data vectors x and y
%           mx, my:     Model vectors x and y. Set to [] for data only.
%           scale:      Scalefactor for arrows (make them bigger or
%                       smaller as desired; use 1 for no scaling)
%           
% OPTIONAL DATA ERROR ELLIPSES
%           Include data uncertainties in order to plot error ellipses. Use
%           either of these forms:
%
%           Option (1): Separate sigma_x, sigma_y (both 1 sigma values), 
%           and optional rotation angle (in radians, CCW from x axis):
%           h = plotarrows(x,y,dx,dy,mx,my,scale,d_sigma_x,d_sigma_y,rot,varargin)
%
%           Option (2): Full covariance matrix (covariance matrix must be
%           in form of E,N,E,N etc.):
%           h = plotarrows(x,y,dx,dy,mx,my,scale,cov,varargin)
%
% OPTIONS FLAGS        
%           confidence:     Default: 2. Plot 1 sigma, 2 sigma, or other
%                           size error ellipses. Note that using a value of
%                           2.55 yields a ~99% confidence ellipse (if
%                           preferred over 68% or 95% etc.). Note, for
%                           Gaussian distributions, use erf(sigma/sqrt(2))
%
%           plotstations:   Default: 1. Plot dots for the stations
%
%           scalebar:       Plot a scale bar of arrow length in the lower-
%                           left corner.
%
%           arrowunits:     Text string with units to use in the scalebar,
%                           above. Does NOT affect any plots, just the
%                           legened text.
%
%           <style options> See flags in code
%

%% Defaults and user input

% Beginning of varargin may be cov, or else three things: sigma_x, sigma_y,
% and rotation. Or none of these!
covtype = 'none';            
k = 1;
if ~isempty(varargin)
    if ~ischar(varargin{1}) % Not a string
        if (size(varargin{1},1) == size(varargin{1},2)) && size(varargin{1},1)>1    %If it's an array
            covtype = 'full';
            cov = varargin{1};
            k = 2;
        else
            covtype = 'sigma';
            sigma_x = varargin{1};
            sigma_y = varargin{2};
            sigma_rot = varargin{3};
            if sigma_rot==0
                sigma_rot = [];
            end
            k = 4;
        end
    end
end

% Anything else in varargin must be user options
options.confidence = 2;             % Plot what sigma ellipses? Default is 2 sigma (95% bounds). Set to 1 for 68% bounds.
options.colors.arrows.model = 'r';
options.colors.arrows.data = 'k';
options.colors.errors = [.6 .6 .6];
options.arrowheads.scale = .2;
options.arrowheads.angle = .33;
options.linewidth = 1;
options.sigma = 1;      % 1 sigma or 2 sigma
confidence_specified = 0;   % For warning text
options.plotstations = 1;
options.colors.stationcolor = 'k';
options.stationsize = 5;
options.scalebar.plot = 0;
options.scalebar.units = [];
% User choices
while k <= length(varargin)
    if ischar(varargin{k})
        switch (lower(varargin{k}))
            case {'confidence'}
                confidence_specified = 1;
                options.confidence = varargin{k+1};
                k = k + 1;
            case {'datacolor'}
                options.colors.arrows.data = varargin{k+1};
                k = k + 1;
            case {'modelcolor'}
                options.colors.arrows.model = varargin{k+1};
                k = k + 1;
            case {'errorcolor'}
                options.colors.errors = varargin{k+1};
                k = k + 1;
            case {'arrowheadscale'}
                options.arrowheads.scale = varargin{k+1};
                k = k + 1;
            case {'arrowheadangle'}
                options.arrowheads.angle = varargin{k+1};
                k = k + 1;
            case {'linewidth'}
                options.linewidth = varargin{k+1};
                k = k + 1;
            case {'plotstations'}
                options.plotstations = varargin{k};
            case {'stationcolor'}
                options.colors.stationcolor = varargin{k+1};
                k = k + 1;
            case {'stationsize'}
                options.stationsize = varargin{k+1};
                k = k + 1;
            case {'scalebar'}
                options.scalebar.plot = varargin{k};
            case {'arrowunits'}
                options.scalebar.units = varargin{k+1};
                k = k + 1;
            otherwise
                warning('Parameter not recognized!')
                disp(varargin{k});
        end
    end
    k = k + 1;
end
% Scale up arrowheads with arrow scale; this isn't perfect but mostly works
options.arrowheads.scale = options.arrowheads.scale * scale;

if ~confidence_specified
    disp('Warning: Plotting default 2 sigma ellipses (~95%)');
end

% Partly from arrows2.m, which is in turn from DM_Quiver.m (Peter's code)

%% Some processing
% Turn everything into row vectors
if size(x,1) > 1; x =x'; end
if size(y,1) > 1; y =y'; end
if size(dx,1) > 1; dx =dx'; end
if size(dy,1) > 1; dy =dy'; end
if size(mx,1) > 1; mx =mx'; end
if size(my,1) > 1; my =my'; end


%% Build arrows
[dX,dY] = arrowfunc(x,y,dx,dy,options.arrowheads.scale,options.arrowheads.angle); % Data vectors
if ~isempty(mx)
    [mX,mY] = arrowfunc(x,y,mx,my,options.arrowheads.scale,options.arrowheads.angle); % Predicted data (model) vectors
else
    mX = [];
    mY = [];
end


%% Error ellipses for data

% Uncertainties are in a chi-squared distribution (not a Gaussian!)
% WAIT! Why would that be true? Wouldn't they just be Gaussian??????
% I think the data errors are simply Gaussian (by definition); only the
% MODEL uncertaintes behave according to the chi-squared distribution. BUT
% note this is NOT the same behavior as we had before using errlps95, which
% uses a chi-squared distribution! I do not really understand this.
% if options.confidence==2    % To avoid using the stats toolbox
%     % confscale = chi2inv(erf(2/sqrt(2)),2);  % Note: erf(2/sqrt(2)) yields 95.45%: The 2 sigma confidence bounds
%     chi2p = 6.1801;
% else
%     P = erf(options.confidence/sqrt(2));
%     chi2p = chi2inv(P,2);
% end

hold on;
if strcmp(covtype,'full')
    % Full covariance matrix was input
    counter = 1;
    for i=1:length(x)        
        % Code taken and modified from errlps95 but NOT USING chi2inv
        [eigenvecs,eigenvals]=eig(full(cov(counter:counter+1,counter:counter+1)));        
        [maxval,index]=max(diag(eigenvals));
        maxvec=eigenvecs(:,index);
        minval=min(diag(eigenvals));
        azmax=atan2(maxvec(2),maxvec(1));
        %         emax=sqrt(chi2p*maxval);      % Chi-squared distribution only!! For MODEL covariances, I think
        %         emin=sqrt(chi2p*minval);
        emax=sqrt(maxval)*options.confidence;       % This yields 1 sigma; convert to other if desired
        emin=sqrt(minval)*options.confidence;       % This yields 1 sigma; convert to other if desired         
        el = ellipse_plot(dX(2,i),dY(2,i), scale*emax, scale*emin, azmax);
        set(el,'Color',options.colors.errors, 'linewidth', options.linewidth);
        counter = counter + 2;
    end

elseif (strcmp(covtype,'sigma'))
    % Here, sigma_x, sigma_y, and possibly a rotation were input
    if isempty(sigma_rot);
        sigma_rot = zeros(size(sigma_x));
    end
    for i = 1:length(x)
        % Scale is arbitrary scaling, the confidence scales from 1 sigma to
        % something else if desired.
        el = ellipse_plot(dX(2,i), dY(2,i), options.confidence*scale*sigma_x(i), options.confidence*scale*sigma_y(i), sigma_rot(i));
        set(el,'Color',options.colors.errors, 'linewidth', options.linewidth);

    end
  
end

%% Plot stations and arrows!
plot(gca, dX,dY, 'color', options.colors.arrows.data, 'linewidth', options.linewidth)
plot(gca, mX,mY, 'color', options.colors.arrows.model, 'linewidth', options.linewidth)
if options.plotstations
    plot(x,y,'linestyle','none','marker','.','color',options.colors.stationcolor,'markersize',options.stationsize)
end

%% Scale bar (optional)
if options.scalebar.plot
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    hh = yl(2)-yl(1);
    ww = xl(2)-xl(1);
    offset_x = .2*ww;
    offset_y = .2*hh;
    
    scalebarsize = options.confidence*scale*round(max(max([dx;dy;mx;my])));
    line([xl(1)+offset_x,(xl(1)+scalebarsize)+offset_x],[yl(1)+offset_y yl(1)+offset_y],'color',options.colors.arrows.data,'linewidth',options.linewidth)
    h = text(xl(1)+offset_x + 1.3*scalebarsize,yl(1)+offset_y,[num2str(scalebarsize), ' ', options.scalebar.units]);
end


%% Helper functions 

    function [X,Y] = arrowfunc(x,y,u,v,alpha,beta)  
        ns = length(x);
        X=[x;x+scale*u;repmat(NaN,1,ns)];
        Y=[y;y+scale*v;repmat(NaN,1,ns)];
        %---Arrow heads---%
        %alpha=0.2; beta=0.33;
        Up=scale*u; Vp=scale*v;
        L=sqrt(sum([X(1,:)-X(2,:);Y(1,:)-Y(2,:)].^2));
        I=find(L>3);
        Up(I)=Up(I)./(L(I)/3);
        Vp(I)=Vp(I)./(L(I)/3);
        X=[X;X(2,:)-alpha*(Up+beta*(Vp+eps));X(2,:);X(2,:)-alpha*(Up-beta*(Vp+eps)); ...
            repmat(NaN,1,ns)];
        Y=[Y;Y(2,:)-alpha*(Vp-beta*(Up+eps));Y(2,:);Y(2,:)-alpha*(Vp+beta*(Up+eps));...
            repmat(NaN,1,ns)];       
    end

end