function stop = outfun(x, optimValues, state)
% this function is called after every iteration of the LFM - it saves the
% current state and also reconstructs the envelopes with the current model
% and prints a .png so we can track progress

% Will Wilkinson - Jan 2018

stop = false;
    global iterResult;
    %{
    x(iterResult.lsig_ind) = sigmoid(x(iterResult.lsig_ind),iterResult.lrng);
    x(iterResult.gsig_ind) = sigmoid(x(iterResult.gsig_ind),iterResult.grng);
    x(iterResult.Dsig_ind) = sigmoid(x(iterResult.Dsig_ind),iterResult.Drng);
    x(iterResult.Ssig_ind) = sigmoid(x(iterResult.Ssig_ind),iterResult.Srng);
    %}
    optimLatest = struct;
    %optimLatest.theta = x;
    optimLatest.theta = iterResult.theta;
    optimLatest.optimValues = optimValues;
    current_dir = pwd;
    if strcmp(current_dir(end-10:end),'experiments')
        save(strcat('../audio_lfm_data/optimLatest_',iterResult.fileName),'-struct','optimLatest');
    else
        save(strcat('audio_lfm_data/optimLatest_',iterResult.fileName),'-struct','optimLatest');
    end
    global iterCount;
    iterCount = optimValues.iteration;
    iterResult.theta_opt = iterResult.theta;
    runSS(iterResult,0,iterResult.ForceFMAP);
    annotation('textbox',[.725 .0 .3 .3],'String',sprintf('Iteration %d',optimValues.iteration),'FitBoxToText','on','FontSize',15,'FontWeight','bold','LineStyle','none')
    if strcmp(current_dir(end-10:end),'experiments')
        print(sprintf('../figures/%s_iter%d',iterResult.fileName,optimValues.iteration),'-dpng')
    else
        print(sprintf('figures/%s_iter%d',iterResult.fileName,optimValues.iteration),'-dpng')
    end
end