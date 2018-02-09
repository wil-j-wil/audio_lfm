
addpath(genpath('../'))

seed = 12345;
rng(seed) % for repeatability

fileName = 'metal'; % name of wav file
suffix = ''; % suffix for .mat file to be saved (in case of multiple versions)
saveResults = 0;    % save the results to a mat file?
matNameLFM = strcat('../audio_lfm_data/lfm_',fileName,'.mat'); % name of the .mat file saved if saveResults = 1
matNameDemod = strcat('../audio_lfm_data/demod_',fileName,'.mat'); % name of demodulation .mat file
load(matNameLFM)

[y,fs] = audioread(strcat('../audio/',fileName,'.wav')); % reads in the file
    
if fs > 24000
    DS = 2; % Down sample the original if desired
else 
    DS = 1;
end
y = y(1:DS:end);
fs = fs/DS;
numFilt =16;

%% Demodulation

demod_struct = load(matNameDemod);
modulators = demod_struct.A;
carriers = demod_struct.C;

metal_x = makeThisThatLong((1:size(modulators,1))/demod_struct.fs,5000);
metal_y = makeThisThatLong(modulators,5000);

fig2=figure(92);clf
[ha, pos] = tight_subplot(1,1,[.1 .1],[.195 .1],[.1 .05]);
set(fig2, 'Units', 'centimeters', 'Position', [6, 20, 11, 4], 'PaperUnits', 'centimeters', 'PaperSize', [11, 4.1]);
%plot((1:size(modulators,1))/demod_struct.fs,modulators)
axes(ha(1))
plot(metal_x,metal_y,'LineWidth',0.275)
title('Amplitude envelopes of a metal impact sound', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 0.3628 0 0.17])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
print(fig2,'-dpdf','demodulation.pdf')


%% LFM
u = lfmData.filterFinal.uAdjusted';
g = log(1+exp(u));
uVar = lfmData.filterFinal.uVar;

newDur = 1000;

metal_x = makeThisThatLong((1:size(modulators,1))/demod_struct.fs,newDur);
metal_y = makeThisThatLong(modulators,newDur);

uVar = makeThisThatLong(uVar,newDur);

u_long = makeThisThatLong(u,newDur);
g_long = makeThisThatLong(g,newDur);

t = 1:size(u_long,1);

color1 = [0 0 1];
color1_ = [0.7 0.7 1];

uncertaintyUp = u_long'+1.96*sqrt(abs(uVar));
uncertaintyLow = u_long'-1.96*sqrt(abs(uVar));
uncertaintyGUp = log(1+exp(uncertaintyUp));
uncertaintyGLow = log(1+exp(uncertaintyLow));

RemoveNoiseFloor = 5; % set to >0 to remove noise floor
[envSynth,~,~] = reconstruct_sig(lfmData,[],RemoveNoiseFloor);
envSynth_long = makeThisThatLong(envSynth,newDur);

fig3=figure(93);clf
[ha, pos] = tight_subplot(1,3,[.01 .02],[.17 .105],[.05 .02]);

%set(fig3, 'Units', 'centimeters', 'Position', [17, 25, 8.5, 12.5], 'PaperUnits', 'centimeters', 'PaperSize', [8.5, 11.75]);
set(fig3, 'Units', 'centimeters', 'Position', [0, 25, 15, 4.75], 'PaperUnits', 'centimeters', 'PaperSize', [15.75, 4.65]);
%subplot(221)
axes(ha(1))
plot(metal_x,metal_y,'LineWidth',0.275)
title('Amplitude envelopes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 0.2328 0 0.17])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
%set(gca,'XTickLabel',[]);

%subplot(222)
axes(ha(2))
%plot(metal_x,g_long,'color',color1)
h1=fill([metal_x fliplr(metal_x)], [uncertaintyGUp fliplr(uncertaintyGLow)], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
hold on
h2=plot(metal_x,g_long,'color',color1,'DisplayName','Mean');
%legend([h1(1) h2(1)])
legend({'Uncertainty','Mean'}, 'interpreter', 'Latex','Position',[0.53 0.752 0.06 0.06],'fontsize',7)
xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'YTickLabel',[]);


title('Predicted latent force', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
axis([0 0.2328 0 2.6])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)

%subplot(223)
axes(ha(3))
plot(metal_x,envSynth_long,'LineWidth',0.275)
title('Reconstructed amplitude envelopes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 0.2328 0 0.17])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'YTickLabel',[]);

%subplot(224)
%axes(ha(2))
%h1=fill([metal_x fliplr(metal_x)], [uncertaintyUp fliplr(uncertaintyLow)], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
%hold on
%h2=plot(metal_x,u_long,'color',color1,'DisplayName','Mean');
%legend({'Uncertainty','Mean'}, 'interpreter', 'Latex','Position',[0.85 0.828 0.06 0.06],'fontsize',6)

%title('Predicted latent force', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
%axis([0 0.2328 -9 4])
%set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
%set(gca,'XTickLabel',[]);

%set(ha(1),'ylabel','Amplitude');

print(fig3,'-dpdf','lfm_metal.pdf')




%%
fileName = 'stim7_applause'; % name of wav file
suffix = ''; % suffix for .mat file to be saved (in case of multiple versions)
saveResults = 0;    % save the results to a mat file?
matNameLFM = strcat('../audio_lfm_data/lfm_',fileName,'.mat'); % name of the .mat file saved if saveResults = 1
matNameDemod = strcat('../audio_lfm_data/demod_',fileName,'.mat'); % name of demodulation .mat file
load(matNameLFM)

[y,fs] = audioread(strcat('../audio/',fileName,'.wav')); % reads in the file
    
if fs > 24000
    DS = 2; % Down sample the original if desired
else 
    DS = 1;
end
y = y(1:DS:end);
fs = fs/DS;
numFilt =16;

addpath(genpath('../gen_model'))
addpath(genpath('../gen_model/GpMat'))
addpath(genpath('../gppad'))
sig_duration = 2;
%%
seed = 121212;%212121;%10101;%121212;%97531;%45678;
rng(seed) % for repeatability
[cascadeModLFM, genLFM] = trainForces(8,sig_duration,{'_5'});


%% LFM
u = lfmData.filterFinal.uAdjusted';
g = log(1+exp(u));
uVar = lfmData.filterFinal.uVar;

newDur = 10000;

highMod = log(1+exp(cascadeModLFM.u+cascadeModLFM.mu));
highMod_long = makeThisThatLong(highMod,newDur);

actEnv = lfmData.demod_struct.A;
app_x = makeThisThatLong((1:size(actEnv,1))/lfmData.demod_struct.fs,newDur);
app_y = makeThisThatLong(actEnv,newDur);

uVar = makeThisThatLong(uVar,newDur);

u_long = makeThisThatLong(u,newDur);
g_long = makeThisThatLong(g,newDur);

lfPrior = [genLFM{1}.fsamp_tr; genLFM{2}.fsamp_tr; genLFM{3}.fsamp_tr]';
lfPrior_long = makeThisThatLong(lfPrior,newDur);

highModPrior = cascadeModLFM.fsamp_tr';
highModPrior_long = makeThisThatLong(highModPrior,newDur);

latentGen = bsxfun(@times,lfPrior,highModPrior);
latentGen_long = makeThisThatLong(latentGen,newDur);

[envSynth,~,sigSynth] = reconstruct_sig(lfmData,latentGen',[],[],sig_duration*lfmData.demod_struct.fs,[],[],zeros(lfmData.demod_struct.D,1));
envSynth_long = makeThisThatLong(envSynth,newDur);
sigSynth_long = makeThisThatLong(sigSynth,newDur);

t = 1:size(u_long,1);

color1 = [0.3 0.3 0.9];
color2 = [0.3 0.9 0.3];
color3 = [0.9 0.3 0.3];
color4 = [0.5 0.5 0.5];
color5 = [0.6 0.3 0.8];
color5 = [0.2 0.2 0.2];
color1_ = [0.7 0.7 1];

fig4=figure(94);clf
[ha, pos] = tight_subplot(3,2,[.07 .05],[.105 .08]);
set(fig4, 'Units', 'centimeters', 'Position', [0, 25, 14, 13.5], 'PaperUnits', 'centimeters', 'PaperSize', [14.4, 12.4]);

x_limit = 2;

axes(ha(1))
plot(app_x,10*app_y,'LineWidth',0.275)
title('Amplitude envelopes', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 x_limit 0 0.26])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(2))
plot(app_x,g_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,g_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,g_long(:,3),'color',color3,'DisplayName','Mean');
title('Predicted latent forces (transformed)', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(3))
p1=plot(app_x,g_long,'color',color4,'DisplayName','Mean');
hold on
p2=plot(app_x,highMod_long,'color',color5,'LineWidth',1);
title('Calculating the high-level modulator', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(4))
plot(app_x,lfPrior_long(:,1),'color',color1,'LineWidth',0.5);
hold on
plot(app_x,lfPrior_long(:,2),'color',color2,'LineWidth',0.5);
plot(app_x,lfPrior_long(:,3),'color',color3,'LineWidth',0.5);
plot(app_x,highModPrior_long,'color',color5,'LineWidth',1.2);
title('Sample from LFM prior and from high-level modulator', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(5))
plot(app_x,latentGen_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,latentGen_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,latentGen_long(:,3),'color',color3,'DisplayName','Mean');
title('Synthesised latent forces', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
xlabel('Time (seconds)', 'interpreter', 'Latex')

axes(ha(6))
plot(app_x,10*envSynth_long,'LineWidth',0.275)
title('Synthesised amplitude envelopes', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 x_limit 0 0.36])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
xlabel('Time (seconds)', 'interpreter', 'Latex')

print(fig4,'-dpdf','gen_model.pdf')


%%
latentGen_edit = 1.25*latentGen_long;

fig5=figure(95);clf
[ha, pos] = tight_subplot(2,2,[.09 .03],[.14 .075],[.09 .06]);
set(fig5, 'Units', 'centimeters', 'Position', [0, 25, 15, 6.5], 'PaperUnits', 'centimeters', 'PaperSize', [15, 6.4]);

font_size = 8;
axes(ha(1))
j_specgram2(lfmData.demod_struct.y,lfmData.demod_struct.fs,0);
title('Spectrogram of applause sound (a)', 'interpreter', 'Latex');
ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
set(gca,'XTickLabel',[]);
axis([0.5 1.5 0 Inf])

axes(ha(3))
plot(app_x,g_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,g_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,g_long(:,3),'color',color3,'DisplayName','Mean');
title('Predicted latent forces (b)', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 1.5 0 2.5])
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size,'XTick',[0.5 1 1.5],'XTickLabel',[0 0.5 1])
xlabel('Time (seconds)', 'interpreter', 'Latex')

axes(ha(2))
j_specgram2(sigSynth,lfmData.demod_struct.fs,0);
title('Spectrogram of synthesised applause (d)', 'interpreter', 'Latex');
%ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
%xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'XTickLabel',[],'XTick',[0.5 1.25 2]);
set(gca,'YTickLabel',[]);
axis([0.5 2 0 Inf])

axes(ha(4))
plot(app_x,latentGen_edit(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,latentGen_edit(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,latentGen_edit(:,3),'color',color3,'DisplayName','Mean');
title('Synthesised latent forces (c)', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 2 0 2.5])
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'YTickLabel',[],'XTick',[0.5 1.25 2],'XTickLabel',[0 0.5 1]);

print(fig5,'-dpdf','gen_model_spec.pdf','-r1000')




%%
latentGen_edit = 1.25*latentGen_long;

fig0=figure(90);clf
[ha, pos] = tight_subplot(2,1,[.09 .03],[.075 .05],[.11 .03]);
set(fig0, 'Units', 'centimeters', 'Position', [0, 25, 15, 15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15]);

font_size = 10;
axes(ha(1))
j_specgram2(lfmData.demod_struct.y,lfmData.demod_struct.fs,0,1/5);
hold on
plot([0.5 1.5 1.5 0.5 0.5],[1900 1900 2300 2300 1900],'color',[1 0.1 0.1],'LineWidth',1)
title('Spectrogram of applause sound', 'interpreter', 'Latex');
ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
set(gca,'XTickLabel',[]);
axis([0.5 1.5 0 Inf])

axes(ha(2))
plot(app_x,app_y(:,[1:9 11:end]),'color',[0.3 0.3 0.3]);
hold on
plot(app_x,app_y(:,10),'color',[1 0.1 0.1],'LineWidth',1);
title('Amplitude envelopes', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 1.5 0 Inf])
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size,'XTick',[0.5 1 1.5],'XTickLabel',[0 0.5 1])
xlabel('Time (seconds)', 'interpreter', 'Latex')



print(fig0,'-dpng','gen_model_web.png','-r1000')



%%
highMod = log(1+exp(cascadeModLFM.u+cascadeModLFM.mu));
highMod_long = makeThisThatLong(highMod,newDur);

lfPrior = [genLFM{1}.fsamp_tr; genLFM{2}.fsamp_tr; genLFM{3}.fsamp_tr]';
lfPrior_long = makeThisThatLong(lfPrior,newDur);

highModPrior = cascadeModLFM.fsamp_tr';
highModPrior_long = makeThisThatLong(highModPrior,newDur);

color1 = [0.2 0.2 0.9];
color2 = [0.15 0.9 0.45];
color3 = [0.9 0.3 0.1];

fig10=figure(100);clf
[ha, pos] = tight_subplot(5,1,[.04 .03],[.081 .042],[.09 .06]);
set(fig10, 'Units', 'centimeters', 'Position', [0, 25, 12, 8], 'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);

font_size = 6;
axes(ha(1))
j_specgram2(lfmData.demod_struct.y,lfmData.demod_struct.fs,0);
title('(a) Spectrogram of applause sound', 'interpreter', 'Latex');
ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
set(gca,'XTick',[0.5 1 1.5],'XTickLabel',[]);
axis([0.5 1.5 0 Inf])

axes(ha(2))
plot(app_x,g_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,g_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,g_long(:,3),'color',color3,'DisplayName','Mean');
plot(app_x,highMod_long,'k-.','LineWidth',0.9)
title('(b) Predicted latent forces (colours) learnt from (a), and the high-level modulator (black)', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 1.5 0 2.5])
set(gca,'XTickLabel',[],'TickLabelInterpreter', 'latex','fontsize',font_size,'XTick',[0.5 1 1.5])

axes(ha(3))
plot(app_x,lfPrior_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,lfPrior_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,lfPrior_long(:,3),'color',color3,'DisplayName','Mean');
plot(app_x,highModPrior_long,'k-.','LineWidth',0.9)
title('(c) Sample from the prior (multiply colours by black to get (d))', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 2 0 2.5])
set(gca,'XTickLabel',[],'TickLabelInterpreter', 'latex','fontsize',font_size,'XTick',[0.5 1.25 2])

axes(ha(4))
plot(app_x,latentGen_edit(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,latentGen_edit(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,latentGen_edit(:,3),'color',color3,'DisplayName','Mean');
title('(d) Synthesised latent forces', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0.5 2 0 2.5])
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
set(gca,'XTick',[0.5 1.25 2],'XTickLabel',[]);

axes(ha(5))
j_specgram2(sigSynth,lfmData.demod_struct.fs,0,1.5);
title('(e) Spectrogram of synthesised applause, generated by passing (d) through the model', 'interpreter', 'Latex');
%ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
%xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'XTick',[0.5 1.25 2],'XTickLabel',[0 0.5 1]);
axis([0.5 2 0 Inf])
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Frequency (Hz)', 'interpreter', 'Latex');

print(fig10,'-dpdf','gen_model_story.pdf','-r1000')


%%
fileName = 'squeaky_stairs'; % name of wav file
suffix = '3'; % suffix for .mat file to be saved (in case of multiple versions)
saveResults = 0;    % save the results to a mat file?
matNameLFM = strcat('../audio_lfm_data/lfm_',fileName,suffix,'.mat'); % name of the .mat file saved if saveResults = 1
matNameDemod = strcat('../audio_lfm_data/demod_',fileName,'.mat'); % name of demodulation .mat file
load(matNameLFM)

[y,fs] = audioread(strcat('../audio/',fileName,'.wav')); % reads in the file
    
if fs > 24000
    DS = 2; % Down sample the original if desired
else 
    DS = 1;
end
y = y(1:DS:end);
fs = fs/DS;
numFilt =16;

addpath(genpath('../gen_model'))
addpath(genpath('../gen_model/GpMat'))
addpath(genpath('../gppad'))
sig_duration = 2;
%%
seed = 12345;%45678;
rng(seed) % for repeatability
[cascadeModLFM, genLFM] = trainForces(7,sig_duration,{'_5'});


%% LFM
u = lfmData.filterFinal.uAdjusted';
g = log(1+exp(u));
uVar = lfmData.filterFinal.uVar;

newDur = 1000;

highMod = log(1+exp(cascadeModLFM.u+cascadeModLFM.mu));
highMod_long = makeThisThatLong(highMod,newDur);

actEnv = lfmData.demod_struct.A;
app_x = makeThisThatLong((1:size(actEnv,1))/lfmData.demod_struct.fs,newDur);
app_y = makeThisThatLong(actEnv,newDur);

uVar = makeThisThatLong(uVar,newDur);

u_long = makeThisThatLong(u,newDur);
g_long = makeThisThatLong(g,newDur);

lfPrior = [genLFM{1}.fsamp_tr; genLFM{2}.fsamp_tr; genLFM{3}.fsamp_tr]';
lfPrior_long = makeThisThatLong(lfPrior,newDur);

highModPrior = cascadeModLFM.fsamp_tr';
highModPrior_long = makeThisThatLong(highModPrior,newDur);

latentGen = bsxfun(@times,lfPrior,highModPrior);
latentGen_long = makeThisThatLong(latentGen,newDur);

[envSynth,~,sigSynth] = reconstruct_sig(lfmData,latentGen',[],[],sig_duration*lfmData.demod_struct.fs,[],[],zeros(lfmData.demod_struct.D,1));
envSynth_long = makeThisThatLong(envSynth,newDur);
sigSynth_long = makeThisThatLong(sigSynth,newDur);

t = 1:size(u_long,1);

color1 = [0.3 0.3 0.9];
color2 = [0.3 0.9 0.3];
color3 = [0.9 0.3 0.3];
color4 = [0.5 0.5 0.5];
color5 = [0.6 0.3 0.8];
color5 = [0.2 0.2 0.2];
color1_ = [0.7 0.7 1];

fig6=figure(96);clf
[ha, pos] = tight_subplot(3,2,[.07 .05],[.105 .08]);
set(fig6, 'Units', 'centimeters', 'Position', [0, 25, 14, 13.5], 'PaperUnits', 'centimeters', 'PaperSize', [14.4, 12.4]);

x_limit = 2;

axes(ha(1))
plot(app_x,10*app_y,'LineWidth',0.275)
title('Amplitude envelopes', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 x_limit 0 0.26])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(2))
plot(app_x,g_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,g_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,g_long(:,3),'color',color3,'DisplayName','Mean');
title('Predicted latent forces (transformed)', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(3))
p1=plot(app_x,g_long,'color',color4,'DisplayName','Mean');
hold on
p2=plot(app_x,highMod_long,'color',color5,'LineWidth',1);
title('Calculating the high-level modulator', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(4))
plot(app_x,lfPrior_long(:,1),'color',color1,'LineWidth',0.5);
hold on
plot(app_x,lfPrior_long(:,2),'color',color2,'LineWidth',0.5);
plot(app_x,lfPrior_long(:,3),'color',color3,'LineWidth',0.5);
plot(app_x,highModPrior_long,'color',color5,'LineWidth',1.2);
title('Sample from LFM prior and from high-level modulator', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
set(gca,'XTickLabel',[]);

axes(ha(5))
plot(app_x,latentGen_long(:,1),'color',color1,'DisplayName','Mean');
hold on
plot(app_x,latentGen_long(:,2),'color',color2,'DisplayName','Mean');
plot(app_x,latentGen_long(:,3),'color',color3,'DisplayName','Mean');
title('Synthesised latent forces', 'interpreter', 'Latex')
axis([0 x_limit 0 3.8])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
xlabel('Time (seconds)', 'interpreter', 'Latex')

axes(ha(6))
plot(app_x,10*envSynth_long,'LineWidth',0.275)
title('Synthesised amplitude envelopes', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 x_limit 0 0.36])
set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
xlabel('Time (seconds)', 'interpreter', 'Latex')

print(fig6,'-dpdf','gen_model_stairs.pdf')

%%

fig7=figure(97);clf
[ha, pos] = tight_subplot(2,1,[.07 .05],[.105 .08]);
set(fig7, 'Units', 'centimeters', 'Position', [0, 25, 19.5, 13], 'PaperUnits', 'centimeters', 'PaperSize', [14.4, 12.4]);

axes(ha(1))
j_specgram2(lfmData.demod_struct.y,lfmData.demod_struct.fs,0);
title(['Orig'],'FontSize',12);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('Time (sec)','FontSize',10);
%set(gca,'YScale','log')


axes(ha(2))
j_specgram2(10*sigSynth,lfmData.demod_struct.fs,0);
title(['Synth'],'FontSize',12);
ylabel('Frequency (Hz)','FontSize',10);
xlabel('Time (sec)','FontSize',10);
%set(gca,'YScale','log')

print(fig7,'-dpdf','gen_model_spec_stairs.pdf')



%% Obj. results

%% List of all examples
fileNames = {'metal','bike_horn2','stim360_cymbal_crash','stim186_hammering','gunshot2','squeaky_stairs2','squeaky_stairs3',...
             'stim7_applause','stim41_breathing','stim162_frog_croaking','glass','dog','stim117_dog_drinking','stim114_dishes_clanking',...
             'stim18_basketball_dribbling','stim373_train_warning_bell','stim132_water_dripping','stim373_train_warning_bell2',...
             'rain','stim399_walking_on_hard_surface','stim78_chimes_in_the_wind','stim211_keys_jingling','stim312_wind',...
             'stim50_camera_snapping_photos','cello','stim398_walking_on_gravel','stim35_boiling_water'}; % name of wav file

objResults = cell(length(fileNames),1);

% Get the objective results and store
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

err = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.errLFM-objResults{k}.errNMF objResults{k}.errNMF-objResults{k}.errNMF objResults{k}.errTNMF-objResults{k}.errNMF] ./ objResults{k}.errNMF)];
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

fig8=figure(98);clf
[ha, pos] = tight_subplot(1,2,[.07 .05],[.21 .11],[0.03 0.02]);
set(fig8, 'Units', 'centimeters', 'Position', [0, 25, 11, 2.7], 'PaperUnits', 'centimeters', 'PaperSize', [11, 2.8]);
axes(ha(1))
%errorbar([1:3],rms_LFM_mean,rms_LFM_down,rms_LFM_up)
plot([0 1 2 3 4],[0 0 0 0 0],'k--','color',[0.35 0.35 0.35],'LineWidth',0.7)
hold on
errorbar([1:3],rms_LFM_median,rms_LFM_down,rms_LFM_up,'x','MarkerSize',1.5,'LineWidth',0.5,'color',[0.1 0 0.8])
errorbar([1:3],rms_tNMF_median,rms_tNMF_down,rms_tNMF_up,'x','MarkerSize',1.5,'LineWidth',0.5,'color',[0.1 0.9 0.35])
axis([0.5 3.5 -1 5.5])
xlabel('Number of latent functions', 'interpreter', 'Latex')
title('RMS error relative to NMF', 'interpreter', 'Latex')
%legend({'NMF','LFM','tNMF'}, 'interpreter', 'Latex','Position',[0.365 0.765 0.04 0.06],'fontsize',4,'linewidth',0.25)
set(gca,'TickLabelInterpreter', 'latex','fontsize',6)
set(gca,'XTick',[1 2 3])

err = cell(4,1);
for k=1:length(fileNames)
    err{objResults{k}.nlf} = [err{objResults{k}.nlf}; ([objResults{k}.lfmdist-objResults{k}.nmfdist objResults{k}.nmfdist-objResults{k}.nmfdist objResults{k}.tnmfdist-objResults{k}.nmfdist] ./ objResults{k}.nmfdist)];
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

axes(ha(2))
%errorbar([1:3],rms_LFM_mean,rms_LFM_down,rms_LFM_up)
plot([0 1 2 3 4],[0 0 0 0 0],'k--','color',[0.35 0.35 0.35],'LineWidth',0.7)
hold on
errorbar([1:3],rms_LFM_median,rms_LFM_down,rms_LFM_up,'x','MarkerSize',1.5,'LineWidth',0.5,'color',[0.1 0 0.8])
errorbar([1:3],rms_tNMF_median,rms_tNMF_down,rms_tNMF_up,'x','MarkerSize',1.5,'LineWidth',0.5,'color',[0.1 0.9 0.35])
axis([0.5 3.5 -1 1])
set(gca,'TickLabelInterpreter', 'latex','fontsize',6)
xlabel('Number of latent functions', 'interpreter', 'Latex')
title('Cosine distance relative to NMF', 'interpreter', 'Latex')
legend({'NMF','LFM','tNMF'}, 'interpreter', 'Latex','Position',[0.596 0.72 0.04 0.06],'fontsize',4,'linewidth',0.25)
%legend({'NMF','LFM','tNMF'}, 'interpreter', 'Latex','Position',[0.608 0.79 0.04 0.06],'fontsize',6)
set(gca,'XTick',[1 2 3])

print(fig8,'-dpdf','objective_results.pdf','-r600')



%%

filelocation = '';
filename = {'Applause','Basketball','Bell2','BikeHorn','Breathing','Camera','Chimes','Cymbal','Dishes','Dog','Frog',...
            'Glass','Gravel','Gunshot','Hammering','Keys','Metal','Stairs2','Stairs3','Wind'};
filename_plot = {'Applause','Basketball','Bell','Bike Horn','Breathing','Camera','Chimes','Cymbal','Dishes','Dog','Frog',...
            'Glass','Gravel','Gunshot','Hammering','Keys','Metal','Stairs','Stairs','Wind'};
nlf = {3, 1, 2, 3, 2, 3, 3, 1, 3, 2, 2,...
        2, 3, 1, 1, 2, 1, 2, 3, 3};
nlf_noncell = [3 1 2 3 2 3 3 1 3 2 2 2 3 1 1 2 1 2 3 3];
startRow = 2;

LFM_ALL = [];
for f_ind=1:length(filename)
    file_to_import = strcat(filelocation,filename{f_ind},'-default-ratings.csv');
    [file_keys,LFM1,LFM2,NMF1,NMF2,TNMF1,TNMF2,anchor] = importfile(file_to_import, startRow);
    LFM_ALL(end+1) = mean([LFM1; LFM2]);
end
[~,I] = sort(LFM_ALL,'descend');
filename = filename(I);
filename_plot = filename_plot(I);
nlf = nlf(I);
nlf_noncell = nlf_noncell(I);

LFM_mean{f_ind} = cell(0);
LFM_med{f_ind} = cell(0);
NMF_mean{f_ind} = cell(0);
NMF_med{f_ind} = cell(0);
tNMF_mean{f_ind} = cell(0);
tNMF_med{f_ind} = cell(0);
fig9=figure(99);clf
[ha, pos] = tight_subplot(1,3,[.07 .012],[.06 .07],[0.07 0.04]);
set(fig9, 'Units', 'centimeters', 'Position', [0, 25, 11, 3.8], 'PaperUnits', 'centimeters', 'PaperSize', [11, 4]);
for f_ind=1:length(filename)
    file_to_import = strcat(filelocation,filename{f_ind},'-default-ratings.csv');
    [file_keys,LFM1,LFM2,NMF1,NMF2,TNMF1,TNMF2,anchor] = importfile(file_to_import, startRow);
    LFM_mean{f_ind} = {mean([LFM1; LFM2]), std([LFM1; LFM2])};
    LFM_med{f_ind} = {median([LFM1; LFM2]), median([LFM1; LFM2])-prctile([LFM1; LFM2],25), prctile([LFM1; LFM2],75)-median([LFM1; LFM2])};
    NMF_mean{f_ind} = {mean([NMF1; NMF2]), std([NMF1; NMF2])};
    NMF_med{f_ind} = {median([NMF1; NMF2]), median([NMF1; NMF2])-prctile([NMF1; NMF2],25), prctile([NMF1; NMF2],75)-median([NMF1; NMF2])};
    tNMF_mean{f_ind} = {mean([TNMF1; TNMF2]), std([TNMF1; TNMF2])};
    tNMF_med{f_ind} = {median([TNMF1; TNMF2]), median([TNMF1; TNMF2])-prctile([TNMF1; TNMF2],25), prctile([TNMF1; TNMF2],75)-median([TNMF1; TNMF2])};
    axes(ha(nlf{f_ind}))
    %subplot(1,3,nlf{f_ind})
    if strcmp(filename{f_ind}, 'Wind')
        pl(f_ind) = plot([LFM_mean{f_ind}{1} NMF_mean{f_ind}{1} tNMF_mean{f_ind}{1}],'ko-','MarkerSize',1.6,'LineWidth',.2);
    else
        pl(f_ind) = plot([LFM_mean{f_ind}{1} NMF_mean{f_ind}{1} tNMF_mean{f_ind}{1}],'o-','MarkerSize',1.6,'LineWidth',.2);
    end
    %errorbar([1:3],[LFM_med{f_ind}{1} NMF_med{f_ind}{1} tNMF_med{f_ind}{1}],[LFM_med{f_ind}{2} NMF_med{f_ind}{2} tNMF_med{f_ind}{2}],[LFM_med{f_ind}{3} NMF_med{f_ind}{3} tNMF_med{f_ind}{3}])
    hold on
    axis([0.5 3.5 0 1])
    set(gca,'XTick',[1 2 3],'XTickLabels',{'LFM','NMF','tNMF'})
    set(gca,'TickLabelInterpreter', 'latex','fontsize',6)
end
axes(ha(1))
ylabel('Realism rating', 'interpreter', 'Latex')
title('1 latent function', 'interpreter', 'Latex')
n_ind = find(nlf_noncell==1);
[~,icons]=legend(pl(n_ind),{filename_plot{n_ind}}, 'interpreter', 'Latex','fontsize',4,'Position',[0.12 0.81 0 0],'color','none');
legend boxoff
numLines = 5;
for i=1:numLines
    LineData = get(icons(2*(i-1)+numLines+1),'XData');
    MarkerData = get(icons(2*(i-1)+numLines+2),'XData');
    NewData6 = [LineData(1)+.2 LineData(2)-.01];
    NewData7 = [MarkerData(1)+.1];
    %// Apply the changes
    set(icons(i),'FontSize',4,'interpreter','latex')
    set(icons(2*(i-1)+numLines+1),'XData',NewData6,'LineWidth',.1)
    set(icons(2*(i-1)+numLines+2),'XData',NewData7,'LineWidth',.1,'MarkerSize',1)
end
axes(ha(2))
title('2 latent functions', 'interpreter', 'Latex')
n_ind = find(nlf_noncell==2);
[~,icons]=legend(pl(n_ind),{filename_plot{n_ind}}, 'interpreter', 'Latex','fontsize',4,'Position',[0.42 0.762 0 0],'color','none');
legend boxoff
numLines = 7;
for i=1:numLines
    LineData = get(icons(2*(i-1)+numLines+1),'XData');
    MarkerData = get(icons(2*(i-1)+numLines+2),'XData');
    NewData6 = [LineData(1)+.2 LineData(2)-.01];
    NewData7 = [MarkerData(1)+.1];
    %// Apply the changes
    set(icons(i),'FontSize',4,'interpreter','latex')
    set(icons(2*(i-1)+numLines+1),'XData',NewData6,'LineWidth',.1)
    set(icons(2*(i-1)+numLines+2),'XData',NewData7,'LineWidth',.1,'MarkerSize',1)
end
set(gca,'YTickLabel',[]);
axes(ha(3))
title('3 latent functions', 'interpreter', 'Latex')
n_ind = find(nlf_noncell==3);
[h3,icons]=legend(pl(n_ind),{filename_plot{n_ind}}, 'interpreter', 'Latex','fontsize',4,'Position',[0.72 0.74 0 0],'color','none');
legend boxoff
numLines = 8;
for i=1:numLines
    LineData = get(icons(2*(i-1)+numLines+1),'XData');
    MarkerData = get(icons(2*(i-1)+numLines+2),'XData');
    NewData6 = [LineData(1)+.2 LineData(2)-.01];
    NewData7 = [MarkerData(1)+.1];
    %// Apply the changes
    set(icons(i),'FontSize',4,'interpreter','latex')
    set(icons(2*(i-1)+numLines+1),'XData',NewData6,'LineWidth',.1)
    set(icons(2*(i-1)+numLines+2),'XData',NewData7,'LineWidth',.1,'MarkerSize',1)
end
set(gca,'YTickLabel',[]);


print(fig9,'-dpdf','listening_test.pdf','-r600')
%print(fig9,'-dpdf','listening_test_dummy.pdf','-r600')



%%
fileName = 'stim78_chimes_in_the_wind'; % name of wav file
suffix = ''; % suffix for .mat file to be saved (in case of multiple versions)
saveResults = 0;    % save the results to a mat file?
matNameLFM = strcat('../audio_lfm_data/lfm_',fileName,'.mat'); % name of the .mat file saved if saveResults = 1
matNameDemod = strcat('../audio_lfm_data/demod_',fileName,'.mat'); % name of demodulation .mat file
load(matNameLFM)

[y,fs] = audioread(strcat('../audio/',fileName,'.wav')); % reads in the file
    
if fs > 24000
    DS = 2; % Down sample the original if desired
else 
    DS = 1;
end
y = y(1:DS:end);
fs = fs/DS;
numFilt =16;

addpath(genpath('../GpSamp'))
addpath(genpath('../GpSamp/GpMat'))
addpath(genpath('../gppad'))
sig_duration = 2;

%%
newDur = 1000;

actEnv = lfmData.demod_struct.A;
app_x = makeThisThatLong((1:size(actEnv,1))/lfmData.demod_struct.fs,newDur);
app_y = makeThisThatLong(actEnv,newDur);
sub_ind = 12;
cen_freq = lfmData.demod_struct.Fc(sub_ind);

fig0=figure(90);clf
[ha, pos] = tight_subplot(2,1,[.08 .03],[.095 .05],[.142 .03]);
set(fig0, 'Units', 'centimeters', 'Position', [0, 25, 15, 15], 'PaperUnits', 'centimeters', 'PaperSize', [15, 15]);

font_size = 12;
axes(ha(1))
j_specgram2(lfmData.demod_struct.y,lfmData.demod_struct.fs,0,1);
hold on
plot([0 2 2 0 0],[cen_freq-200 cen_freq-200 cen_freq+200 cen_freq+200 cen_freq-200],'color',[1 0.1 0.1],'LineWidth',1.5)
title('Spectrogram of wind chimes', 'interpreter', 'Latex');
ylabel('Frequency (Hz)', 'interpreter', 'Latex');
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size)
set(gca,'XTick',[0 0.5 1 1.5 2],'XTickLabel',[]);
axis([0 2 0 Inf])

axes(ha(2))
plot(app_x,app_y(:,[1:9 11:end]),'color',[0.3 0.3 0.3]);
hold on
plot(app_x,app_y(:,sub_ind),'color',[1 0.1 0.1],'LineWidth',1.5);
title('Amplitude envelopes', 'interpreter', 'Latex')
%axis([0 x_limit 0 3.5])
axis([0 2 0 Inf])
set(gca,'TickLabelInterpreter', 'latex','fontsize',font_size,'XTick',[0 0.5 1 1.5 2],'XTickLabel',[0 0.5 1 1.5 2])
xlabel('Time (seconds)', 'interpreter', 'Latex')



print(fig0,'-dpng','gen_model_web.png','-r100')



%% LFM
[w,h] = nnmf(metal_y,3);
%%
modulators = lfmData.demod_struct.A;

u = lfmData.filterFinal.uAdjusted';
g = log(1+exp(u));
uVar = lfmData.filterFinal.uVar;

newDur = 1000;

metal_x = makeThisThatLong((1:size(modulators,1))/demod_struct.fs,newDur);
metal_y = makeThisThatLong(modulators,newDur);


uVar = makeThisThatLong(uVar,newDur);

u_long = makeThisThatLong(u,newDur);
g_long = makeThisThatLong(g,newDur);

t = 1:size(u_long,1);

color1 = [0 0 1];
color2 = [1 0 0];
color3 = [0 1 0];
color1_ = [0.7 0.7 1];
color2_ = [1 0.7 0.7];
color3_ = [0.7 1 0.7];

uncertaintyUp = u_long'+1.96*sqrt(abs(uVar));
uncertaintyLow = u_long'-1.96*sqrt(abs(uVar));
uncertaintyGUp = log(1+exp(uncertaintyUp));
uncertaintyGLow = log(1+exp(uncertaintyLow));

RemoveNoiseFloor = 5; % set to >0 to remove noise floor
[envSynth,~,~] = reconstruct_sig(lfmData,[],RemoveNoiseFloor);
envSynth_long = makeThisThatLong(envSynth,newDur);

fig3=figure(93);clf
[ha, pos] = tight_subplot(1,2,[.01 .02],[.17 .105],[.05 .02]);

%set(fig3, 'Units', 'centimeters', 'Position', [17, 25, 8.5, 12.5], 'PaperUnits', 'centimeters', 'PaperSize', [8.5, 11.75]);
set(fig3, 'Units', 'centimeters', 'Position', [0, 25, 30, 9.5], 'PaperUnits', 'centimeters', 'PaperSize', [30, 9.5]);
%subplot(221)
%{
axes(ha(1))
plot(metal_x,metal_y,'LineWidth',0.275)
title('Amplitude envelopes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 2 0 0.03])
set(gca,'TickLabelInterpreter', 'latex','fontsize',14)
%set(gca,'XTickLabel',[]);
%}

%subplot(222)
axes(ha(1))
%plot(metal_x,g_long,'color',color1)
h3=fill([metal_x fliplr(metal_x)], [uncertaintyGUp(3,:) fliplr(uncertaintyGLow(3,:))], color3_, 'edgecolor',color3_,'DisplayName','Uncertainty');
hold on
h1=fill([metal_x fliplr(metal_x)], [uncertaintyGUp(1,:) fliplr(uncertaintyGLow(1,:))], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
h2=fill([metal_x fliplr(metal_x)], [uncertaintyGUp(2,:) fliplr(uncertaintyGLow(2,:))], color2_, 'edgecolor',color2_,'DisplayName','Uncertainty');
plot(metal_x,g_long(:,1),'color',color1,'LineWidth',1.5);
%hold on
plot(metal_x,g_long(:,2),'color',color2,'LineWidth',1.5);
plot(metal_x,g_long(:,3),'color',color3,'LineWidth',1.5);
%legend([h1(1) h2(1)])
%legend({'Uncertainty','Mean'}, 'interpreter', 'Latex','Position',[0.53 0.752 0.06 0.06],'fontsize',7)
xlabel('Time (seconds)', 'interpreter', 'Latex')
set(gca,'YTick',[],'YTickLabel',[]);


title('Predicted latent forces (our model)', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
axis([0 2 0 2.2])
set(gca,'TickLabelInterpreter', 'latex','fontsize',12)

%subplot(223)
%{
axes(ha(3))
plot(metal_x,envSynth_long,'LineWidth',0.275)
title('Reconstructed amplitude envelopes', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 2 0 0.03])
set(gca,'TickLabelInterpreter', 'latex','fontsize',14)
set(gca,'YTickLabel',[]);
%}
axes(ha(2))
plot(metal_x,w(:,1),'color',color1,'LineWidth',1.5)
hold on
plot(metal_x,w(:,2),'color',color2,'LineWidth',1.5)
plot(metal_x,w(:,3),'color',color3,'LineWidth',1.5)
title('NMF Activations', 'interpreter', 'Latex')
xlabel('Time (seconds)', 'interpreter', 'Latex')
%ylabel('Amplitude', 'interpreter', 'Latex')
axis([0 2 0 0.035])
set(gca,'TickLabelInterpreter', 'latex','fontsize',12)
set(gca,'YTick',[],'YTickLabel',[]);

%subplot(224)
%axes(ha(2))
%h1=fill([metal_x fliplr(metal_x)], [uncertaintyUp fliplr(uncertaintyLow)], color1_, 'edgecolor',color1_,'DisplayName','Uncertainty');
%hold on
%h2=plot(metal_x,u_long,'color',color1,'DisplayName','Mean');
%legend({'Uncertainty','Mean'}, 'interpreter', 'Latex','Position',[0.85 0.828 0.06 0.06],'fontsize',6)

%title('Predicted latent force', 'interpreter', 'Latex')
%xlabel('Time (seconds)', 'interpreter', 'Latex')
%axis([0 0.2328 -9 4])
%set(gca,'TickLabelInterpreter', 'latex','fontsize',8)
%set(gca,'XTickLabel',[]);

%set(ha(1),'ylabel','Amplitude');

%print(fig3,'-dpdf','lfm_metal.pdf')
print(fig3,'-dpng','lfm_vs_nmf_web.png','-r100')