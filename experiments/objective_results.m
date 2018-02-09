%% Calculate and plot the objective results (Euclidean and cosine distance of reconstruction to original signal) of LFM vs. NMF vs. tNMF

%% List of all examples
fileNames = {'metal','bike_horn2','stim360_cymbal_crash','stim186_hammering','gunshot2','squeaky_stairs2','squeaky_stairs3',...
             'stim7_applause','stim41_breathing','stim162_frog_croaking','glass','dog','stim117_dog_drinking','stim114_dishes_clanking',...
             'stim18_basketball_dribbling','stim373_train_warning_bell','stim132_water_dripping','stim373_train_warning_bell2',...
             'rain','stim399_walking_on_hard_surface','stim78_chimes_in_the_wind','stim211_keys_jingling','stim312_wind',...
             'stim50_camera_snapping_photos','cello','stim398_walking_on_gravel','stim35_boiling_water'}; % name of wav file

objResults = cell(length(fileNames),1);

%% Get the objective results and store
for i=1:length(fileNames)
    fprintf('retrieving objective results: %s\n',fileNames{i})
    load(strcat('../audio_lfm_data/lfm_',fileNames{i},'.mat'));
    objResults{i} = lfmData.objective_results;
    objResults{i}.fileName = fileNames{i};
    if strcmp(objResults{i}.fileName(1:min(4,end)),'stim')
        undersc = strfind(objResults{i}.fileName,'_');
        objResults{i}.fileName = objResults{i}.fileName(undersc+1:end);
    end
    objResults{i}.nlf = lfmData.lfmOpt.nlf;
end

%%
figure(7);clf
for j=1:size(objResults,1)
    if objResults{j}.errLFM >= objResults{j}.errNMF;
        err = 100*(objResults{j}.errLFM/objResults{j}.errNMF - 1);
    else
        err = 100*(1 - objResults{j}.errNMF/objResults{j}.errLFM);
    end
    if objResults{j}.lfmdist >= objResults{j}.nmfdist;
        cdist = 100*(objResults{j}.lfmdist/objResults{j}.nmfdist - 1);
    else
        cdist = 100*(1 - objResults{j}.nmfdist/objResults{j}.lfmdist);
    end
    if objResults{j}.nlf == 1
        p1=plot(err,cdist,'b.','markersize',20);
    elseif objResults{j}.nlf == 2
        p2=plot(err,cdist,'g.','markersize',20);
    elseif objResults{j}.nlf == 3
        p3=plot(err,cdist,'r.','markersize',20);
    else
        p4=plot(err,cdist,'c.','markersize',20);
    end
    hold on
    text(err-5,cdist-5,objResults{j}.fileName,'interpreter','none','FontSize',11)
end
legend([p1 p2 p3 p4],{'1 Force','2 Forces','3 Forces','4 Forces'},'Location','northwest')
xlim([-250 250])
ylim([-250 250])
%zeroAxes(gca)
set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
set(gca,'FontSize',12)
%axis off
line([0 0], [250 -250])
line([250 -250], [0 0])
box off
xlabel('Euclidean Distance (%)')
ylabel('Cosine Distance (%)')
title('LFM vs NMF error comparison')

%% RMS Error bar chart
err = cell(4,1);
fnames = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.errLFM objResults{k}.errNMF objResults{k}.errTNMF] ./ objResults{k}.errNMF)];
    fnames{objResults{k}.nlf}{end+1} = objResults{k}.fileName;
end
map = [0.2, 0.2, 0.9; 0.2, 0.9, 0.2; 0.9, 0.2, 0.2];
figure(8);clf
annotation('textbox',[0.44 0.885 0.1 0.1],'FontSize',16,'FontWeight','Bold','String','RMS Error / Euclidean Distance','FitBoxToText','on', 'EdgeColor', 'none','Color','k')
subplot(221)
bar(err{1})
colormap(map)
ylim([0 2])
title('1 Latent Function')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{1},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(222)
bar(err{2})
ylim([0 2])
title('2 Latent Functions')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{2},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(223)
bar(err{3})
ylim([0 2])
title('3 Latent Functions')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{3},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(2,6,10)
boxplot(err{1}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 7])
title('1 Latent Function')
set(gca, 'XTickLabel', {'LFM','tNMF'})
subplot(2,6,11)
boxplot(err{2}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 7])
title('2 Latent Functions')
set(gca, 'XTickLabel', {'LFM','tNMF'})
subplot(2,6,12)
boxplot(err{3}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 7])
title('3 Latent Functions')
set(gca, 'XTickLabel', {'LFM','tNMF'})
%{
subplot(224)
bar(err{4})
ylim([0 2])
title('4 Latent Functions')
%}
%% Cosine Distance bar chart
err = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.lfmdist objResults{k}.nmfdist objResults{k}.tnmfdist] ./ objResults{k}.nmfdist)];
end
map = [0.2, 0.2, 0.9; 0.2, 0.9, 0.2; 0.9, 0.2, 0.2];
figure(9);clf
annotation('textbox',[0.47 0.885 0.1 0.1],'FontSize',16,'FontWeight','Bold','String','Cosine Distance','FitBoxToText','on', 'EdgeColor', 'none','Color','k')
subplot(221)
bar(err{1})
colormap(map)
ylim([0 2])
title('1 Latent Function')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{1},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(222)
bar(err{2})
ylim([0 2])
title('2 Latent Functions')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{2},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(223)
bar(err{3})
ylim([0 2])
title('3 Latent Functions')
legend('LFM','NMF','tNMF')
set(gca, 'XTickLabel', fnames{3},'XTickLabelRotation',15,'TickLabelInterpreter','none')
subplot(2,6,10)
boxplot(err{1}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 3])
title('1 Latent Function')
set(gca, 'XTickLabel', {'LFM','tNMF'})
subplot(2,6,11)
boxplot(err{2}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 3])
title('2 Latent Functions')
set(gca, 'XTickLabel', {'LFM','tNMF'})
subplot(2,6,12)
boxplot(err{3}(:,[1 3]))
hold on
plot([0 1 2 3],[1 1 1 1],'g','LineWidth',2)
legend('NMF','Location','Northwest')
ylim([0 3])
title('3 Latent Functions')
set(gca, 'XTickLabel', {'LFM','tNMF'})
%{
subplot(224)
bar(err{4})
ylim([0 2])
title('4 Latent Functions')
%}


%%
err = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.errLFM objResults{k}.errNMF objResults{k}.errTNMF] ./ objResults{k}.errNMF)];
end

%rms_LFM_mean = [mean(err{1}(:,1)) mean(err{2}(:,1)) mean(err{3}(:,1))];
rms_LFM_median = [median(err{1}(:,1)) median(err{2}(:,1)) median(err{3}(:,1))];
rms_LFM_up = [prctile(err{1}(:,1),75) prctile(err{2}(:,1),75) prctile(err{3}(:,1),75)];
rms_LFM_up = rms_LFM_up - rms_LFM_median;
rms_LFM_down = [prctile(err{1}(:,1),25) prctile(err{2}(:,1),25) prctile(err{3}(:,1),25)];
rms_LFM_down = rms_LFM_median - rms_LFM_down;
rms_tNMF_median = [median(err{1}(:,3)) median(err{2}(:,3)) median(err{3}(:,3))];
rms_tNMF_up = [prctile(err{1}(:,3),75) prctile(err{2}(:,3),75) prctile(err{3}(:,3),75)];
rms_tNMF_up = rms_tNMF_up - rms_tNMF_median;
rms_tNMF_down = [prctile(err{1}(:,3),25) prctile(err{2}(:,3),25) prctile(err{3}(:,3),25)];
rms_tNMF_down = rms_tNMF_median - rms_tNMF_down;

figure(10);clf
subplot(121)
%errorbar([1:3],rms_LFM_mean,rms_LFM_down,rms_LFM_up)
plot([0 1 2 3 4],[1 1 1 1 1],'k--','LineWidth',2,'color',[0.35 0.35 0.35],'LineWidth',3)
hold on
errorbar([1:3],rms_LFM_median,rms_LFM_down,rms_LFM_up,'o','MarkerSize',10,'LineWidth',3)
errorbar([1:3],rms_tNMF_median,rms_tNMF_down,rms_tNMF_up,'o','MarkerSize',10,'LineWidth',3)
axis([0.5 3.5 0 6.5])

err = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.lfmdist objResults{k}.nmfdist objResults{k}.tnmfdist] ./ objResults{k}.nmfdist)];
end

%rms_LFM_mean = [mean(err{1}(:,1)) mean(err{2}(:,1)) mean(err{3}(:,1))];
rms_LFM_median = [median(err{1}(:,1)) median(err{2}(:,1)) median(err{3}(:,1))];
rms_LFM_up = [prctile(err{1}(:,1),75) prctile(err{2}(:,1),75) prctile(err{3}(:,1),75)];
rms_LFM_up = rms_LFM_up - rms_LFM_median;
rms_LFM_down = [prctile(err{1}(:,1),25) prctile(err{2}(:,1),25) prctile(err{3}(:,1),25)];
rms_LFM_down = rms_LFM_median - rms_LFM_down;
rms_tNMF_median = [median(err{1}(:,3)) median(err{2}(:,3)) median(err{3}(:,3))];
rms_tNMF_up = [prctile(err{1}(:,3),75) prctile(err{2}(:,3),75) prctile(err{3}(:,3),75)];
rms_tNMF_up = rms_tNMF_up - rms_tNMF_median;
rms_tNMF_down = [prctile(err{1}(:,3),25) prctile(err{2}(:,3),25) prctile(err{3}(:,3),25)];
rms_tNMF_down = rms_tNMF_median - rms_tNMF_down;

subplot(122)
%errorbar([1:3],rms_LFM_mean,rms_LFM_down,rms_LFM_up)
plot([0 1 2 3 4],[1 1 1 1 1],'k--','LineWidth',2,'color',[0.35 0.35 0.35],'LineWidth',3)
hold on
errorbar([1:3],rms_LFM_median,rms_LFM_down,rms_LFM_up,'o','MarkerSize',10,'LineWidth',3)
errorbar([1:3],rms_tNMF_median,rms_tNMF_down,rms_tNMF_up,'o','MarkerSize',10,'LineWidth',3)
axis([0.5 3.5 0 2])