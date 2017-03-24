%scna through

%% Changing sp\_Lgauss to Raco\_gauss
cutsize =4;
newprefix = 'AdaRBFGS';
for i=1:length(files)
    pathname = [ path files(i).name ]
    eval(['load ' pathname ]);
    if(length(plotdata.legend) ==4)
        %% Fix labels!
        pause;
        ind1 = 3;
        ind2 = 4;
        leg = plotdata.legend{ind1};
        plotdata.legend{ind1} = [ newprefix  leg(cutsize+1:end) ] ;
        leg = plotdata.legend{ind2};
        plotdata.legend{ind2} = [ newprefix  leg(cutsize+1:end) ] ;
        %% Swap first and last!
        plotdata.legend = [plotdata.legend ;   plotdata.legend{2}];
        plotdata.legend(2) =[];
        plotdata.TIMES = [plotdata.TIMES ;   plotdata.TIMES{2}];
        plotdata.TIMES(2) =[];
        plotdata.ERRORS = [plotdata.ERRORS ;   plotdata.ERRORS{2}];
        plotdata.ERRORS(2) =[];
        plotdata.EXTRAFIELD = [plotdata.EXTRAFIELD ;   plotdata.EXTRAFIELD{2}];
        plotdata.EXTRAFIELD(2) =[];
        %% Correct ShultzNewton to Newton Shulz
        plotdata.legend{end} = 'NewtonSchulz';
        close all;
        opts.ylimits = 'adapt';
        prettyPlot_plotdata(plotdata,opts);
        %  eval(['save(''' pathname ''',''plotdata'');' ])
    end
end
