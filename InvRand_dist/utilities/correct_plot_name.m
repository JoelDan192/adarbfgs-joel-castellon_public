%% correct_plot_name
path = 'plotdata/'; 
cd plotdata/;
files = dir( '*.mat');
cd ..
%% Print all names

for i=1:length(files)
   fprintf(' %2.0d :  %s \n',i, files(i).name);
end
%% Changing sp\_Lgauss to Raco\_gauss
cutsize =4;
newprefix = 'AdaRBFGS';
for i=29:length(files)
     pathname = [ path files(i).name ]
     eval(['load ' pathname ]);
     ind1 = 3;
     ind2 = 4;
     leg = plotdata.legend{ind1};
     plotdata.legend{ind1} = [ newprefix  leg(cutsize+1:end) ] ;
     leg = plotdata.legend{ind2};
     plotdata.legend{ind2} = [ newprefix  leg(cutsize+1:end) ] ;
    close all;
    opts.ylimits = 'adapt';
    prettyPlot_plotdata(plotdata,opts);
  %  eval(['save(''' pathname ''',''plotdata'');' ])
end

