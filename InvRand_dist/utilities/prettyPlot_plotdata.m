function  prettyPlot_plotdata(plotdata,opts)
lO = length(plotdata.ERRORS);
for i =1:lO
    lent(i) = length(plotdata.ERRORS{i});
    lowERROR(i) = plotdata.ERRORS{i}(end);
    flops{i} = plotdata.EXTRAFIELD{i}*(1:1:lent(i)); 
end
markers={'o','+','^','*','s','d','v','+','<','>','x','h','.',...
'^','*','o','p','+','<','h','.','>','x','s','d','v',...
'o','p','+','*','s','d','v','+','<','>','x','h','.'};
options.markersize =8;
options.linewidth =1.5;
options.markers= markers(1:lO);
options.colors = lines(lO);
options.xlabel = 'time (s)';
%options.xlabel = 'iterations';
options.ylabel ='Tr(AX_k)';


if(isfield(opts,'ylimits'))
    options.ylimits = [max(10^(-3),min(lowERROR)*(1-0.25)) inf];
   % options.ylimits = [5*10^(-3) inf];
end
options.logScale = 2;
%options.ylabel ='|x -x*|_{A^TA} / |x*|_{A^TA}';
markerspace = [ceil(lent/7)' ceil(rand(lO,1).*(lent./7)')];
options.markerSpacing = markerspace;
lineStyles= {':','--','-',':','--','-',':','--','-',':','--','-'};
options.lineStyles = lineStyles(1:lO);
% options.legend = legendStr;
% Naming
path = [ 'figures/' plotdata.title ]
FigHandle = figure('Position', [0, 0, 750, 300]);
subplot(1,2,1);  %set(gca, 'Position', [0 0 200 400])
prettyPlot(plotdata.TIMES,plotdata.ERRORS, options);
set(gca,'defaultTextFontName', 'Arial')
set(gca,'fontsize',15)
xlhand = get(gca,'xlabel');  set(xlhand,'fontsize',15) ;
xlhand = get(gca,'ylabel');  set(xlhand,'fontsize',15) ;
ch = get(gcf,'children'); set(ch(1), 'fontsize',12);
%Repeat for flops
subplot(1,2,2);   
options= rmfield(options,'ylabel');
% options= rmfield(options,'legend');
options.legend = plotdata.legend;
options.xlabel = 'flops';
prettyPlot(flops ,plotdata.ERRORS, options);
set(gca,'defaultTextFontName', 'Arial')
set(gca,'fontsize',15)
xlhand = get(gca,'xlabel');  set(xlhand,'fontsize',15) ;
xlhand = get(gca,'ylabel');  set(xlhand,'fontsize',15) ;
ch = get(gcf,'children'); set(ch(1), 'fontsize',12);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 12)
set(gca,'ytick',[])
if (isfield(options,'logScale'))
    eval(['print -dpdf ' path '.pdf' ]);
else
    eval(['print -dpdf ' path '-not-log' '.pdf' ]);
end
end