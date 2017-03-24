    close all;
    figure('visible','off');
    if(isnumeric(options.M0)) options.M0 = num2str(options.M0); end
    Prob.title = [ Prob.title '_M0_' options.M0  ];
    plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
    prettyPlot_plotdata(plotdata, options)