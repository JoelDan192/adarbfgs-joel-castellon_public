function  plot_outputs(OUTPUTs,plotoptions)
%plot_outputs plots the time X relative error in log scale of each output, where each
%cell of OUTPUTs is a structure that contains the output of a method.
%
%
% INPUTS
%-------
% a cell array OUTPUTs where each contains
%   OUTPUTs{i}.times    the time recorded at the end of each iteration
%   OUTPUTs{i}.error    the relative error recorded at the end of each iteration
% Copyright (c) 2014.  Robert Gower.

lO = length(OUTPUTs);
MarkerEdgeColors=lines(lO);  
Markers=['o','+','^','*','s','d','v','+','<','>','x','h','.',...
'^','*','o','p','+','<','h','.','>','x','s','d','v',...
'o','p','+','*','s','d','v','+','<','>','x','h','.'];
leg ={};
h =  figure('visible','on');
set(gca,'defaultTextFontName', 'Arial')
set(gca,'FontSize',20)

for i = 1: lO
    output = OUTPUTs{i};
    if(~isfield(output,'times') || ~isfield(output,'errors') )
        display('No times or errors field in output. Please run QuNICsolve again with opts.plotting= 1');
        return;
    end    
    lt = length(output.times);
    if(i==1)
        mark_spacing =1:1:lt;
    elseif(i==5)
        mark_spacing =1:1:lt;
    else
        mark_spacing = 1:1:lt;
    end  
    % plot marker semilogy
    plot(output.times(mark_spacing), (output.errors(mark_spacing)), [':' Markers(i)],'Color',MarkerEdgeColors(i,:), ... 
   'LineWidth',1.1,'MarkerSize',6, 'MarkerEdgeColor', MarkerEdgeColors(lO-i+1,:));
    hold on;
    leg = [leg , {output.name}];
end

legend(leg,'Interpreter', 'latex')
title(output.opts.problem_title,'Interpreter', 'latex')
xlabel('time (s)','Interpreter', 'latex')
ylabel('error','Interpreter', 'latex')
title_name = OUTPUTs{1}.opts.problem_title;
title_name(ismember(title_name,' ,.:;!/\')) = [];
eval([' print  -djpeg ' title_name '.jpg' ]);
print([title_name '.fig']);
hold off;
