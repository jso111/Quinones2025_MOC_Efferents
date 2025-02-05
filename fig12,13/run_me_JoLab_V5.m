
%% this code plots all the figures for ABR and LFP data from Matt McGinley lab
% I used MATLAB 2023; did not test for other versions


a=1;
%cd('Z:\Users\hsrivastava\Data4JO')
% extract LFP data
load("ABR_and_LFP_Data_HS.mat")
[ALL_Binned_data SessionWise_data ALL_units Pupil_and_Response ALL_Binned_PupilData BFSessionWise]=extract_data_from_LFPS(ABR_and_LFP_Data_HS);

% cal R-square
% find_figure('Rsquare with LOOCV')
% subplot(3,2,4)
find_figure('1d');clf
hold on
[RsqvN1]=getRsquare_helper(ALL_Binned_data(:,1:10,1),-1,ALL_Binned_PupilData);
[RsqvP1]=getRsquare_helper(ALL_Binned_data(:,1:10,2),1,ALL_Binned_PupilData);
[allABR]=getR2_withLOOCV_forABT(ABR_and_LFP_Data_HS.ABR);
[ABRRatio]=getR2_withLOOCV_forABR_peakVI_ratio(ABR_and_LFP_Data_HS.ABR);
xticks([1 2 3 4 5 6 8 9])
xticklabels({'ABR-I','ABR-II','ABR-III','ABR-IV','ABR-V','log(V-I ratio)','LFP-P1','LFP-N1'})
ylabel('cross-validated R2')
set(gca,'TickDir','out')
set(gca, 'box', 'off')
set(gcf,'Position',[500   300   900   750])




% % plot mean LFPs on subplots
 subP=1;
% find_figure ('Means')
% subplot(2,4,8)
% plot_mean_LFPs(-ALL_Binned_data(:,1:10,1),subP,-1,ALL_Binned_PupilData)
% subplot(2,4,4)
% plot_mean_LFPs(ALL_Binned_data(:,1:10,2),subP,1,ALL_Binned_PupilData)
% plot_ABR_amps(ABR_and_LFP_Data_HS.ABR,subP)
find_figure('Supp1b')
 plot_p5p1_ratio(ABR_and_LFP_Data_HS.ABR,subP)
 xlim([.05 .55])
 xticks([.1:.1:.5])
 xticklabels([10:10:50])
 ylim([-2.7 1.2])
  xlabel('% of max pupil size')
  set(gcf,'Position',[500   300   900   750])
  ax=gca;
ax.XAxis.FontSize = 18;
   ax.YAxis.FontSize = 18;
 





find_figure('1b')
subP=0;

plot_mean_LFPs(-ALL_Binned_data(:,1:10,1),subP,-1,ALL_Binned_PupilData,SessionWise_data)
plot_mean_LFPs(ALL_Binned_data(:,1:10,2),subP,1,ALL_Binned_PupilData,SessionWise_data)
plot_ABR_amps(ABR_and_LFP_Data_HS.ABR,subP)
legend({'','LFP N1','','','LFP P1','ABR-I','','ABR-II','','ABR-III','','ABR-IV','','ABR-V','',''},'Location','southeast','Box','off')
set(gca,'TickDir','out')
set(gca, 'box', 'off')
set(gcf,'Position',[500   300   900   750])
xlim([.05 .55])
xticks([.1:.1:.5])
xticklabels([10:10:50])
xlabel('% of max pupil size')
ax=gca;
ax.XAxis.FontSize = 18;
   ax.YAxis.FontSize = 18;

% Plot R-square after local averaging
find_figure('Supp1a')
% find_figure('R2_local_averaging')
[R2_locally_averagedN1] = get_local_averaged_R2(Pupil_and_Response,-1);
[R2_locally_averagedP1] = get_local_averaged_R2(Pupil_and_Response,1);
[R2_locally_averagedABR]=get_R2_ABR_local_averaging(ABR_and_LFP_Data_HS.ABR);
legend({'','LFP N1','','LFP P1','','ABR-I','','ABR-II','','ABR-III','','ABR-IV','','ABR-V'},'Location','northwest','Box','off')
xlim([0 6])
ylim([-.01 .06])
xlabel('trials per bin (#)')
ylabel('R^2 for quadratic fit')
set(gca,'TickDir','out')
set(gca, 'box', 'off')
set(gcf,'Position',[500   300   900   750])
ax=gca;
ax.XAxis.FontSize = 18;
   ax.YAxis.FontSize = 18;



% plot LFP population time series
% subplot(3,2,3)
find_figure('1c');clf
% find_figure('LFP TimeSeries')
plot_lfp_timeseries_pupil_binned(ABR_and_LFP_Data_HS.LFPtimseries)
set(gcf,'Position',[500   300   900   750])
   

% plot ABR  time series pop data
% subplot(3,2,1)
find_figure('1a');clf
% find_figure('ABR TimeSeries')
plot_ABR_timeseries(ABR_and_LFP_Data_HS.ABRtimeseries)
set(gcf,'Position',[500   300   900   750])




%% all helper function below


function[ALL_Binned_data SessionWise_data ALL_units Pupil_and_Response ALL_Binned_PupilData BFSessionWise]=extract_data_from_LFPS(ABR_and_LFP_Data_HS);

inclusion_criteria=-1; %-1 for based on LFPs, 0 if spikes
unit_cnt=1;
unit_session=1;
sessionCNT=1;
for fols=1:length(ABR_and_LFP_Data_HS.LFP)
    unit_session=1;
    % ALL_data=[];ALL_units=[];

    resp_channels=[];resp_channels=unique(ABR_and_LFP_Data_HS.LFP(fols).BFfromSpikes(:,2));


    upper_boundary=[];upper_boundary=ABR_and_LFP_Data_HS.LFP(fols).ClusterID(find(ABR_and_LFP_Data_HS.LFP(fols).ClusterID(:,3)<=ABR_and_LFP_Data_HS.LFP(fols).IC_span(2),1,"last"),1);



    CN_DC_boundary=[];CN_DC_boundary=upper_boundary(1)-300; % in microns

    channelNumAtCN_DC_boundary=[];channelNumAtCN_DC_boundary=ABR_and_LFP_Data_HS.LFP(fols).ClusterID(find(ABR_and_LFP_Data_HS.LFP(fols).ClusterID(:,1)<CN_DC_boundary,1,"last"),3);



    for num_sites=1:size(ABR_and_LFP_Data_HS.LFP(fols).PeakN1,1)


        if sum(ismember(resp_channels,ABR_and_LFP_Data_HS.LFP(fols).ChannelID(num_sites)))>inclusion_criteria % CHANGE TO 0 FOR INCLUSION CRITRIA BASED ON SPIKES
            if ABR_and_LFP_Data_HS.LFP(fols).ChannelID(num_sites)<=channelNumAtCN_DC_boundary

                PupilBinEdges=[];PupilBinEdges=ABR_and_LFP_Data_HS.LFP(fols).PupilBinEdges;
                PupilSizeFreqXReps=[];PupilSizeFreqXReps=ABR_and_LFP_Data_HS.LFP(fols).PupilSizeFreqXReps;
                WalkBins=[];WalkBins=ABR_and_LFP_Data_HS.LFP(fols).WalkBins;
                



                current_siteN1=[];current_siteP1=[];
                current_siteN1=squeeze(ABR_and_LFP_Data_HS.LFP(fols).PeakN1(num_sites,:,:));
                current_siteP1=squeeze(ABR_and_LFP_Data_HS.LFP(fols).PeakP1(num_sites,:,:));
                non_responding_freqsP1=find(isnan(current_siteP1(:,1))==1);

                non_responding_freqsN1=find(isnan(current_siteN1(:,1))==1);
                if length(non_responding_freqsN1)<26 & length(non_responding_freqsP1)<26
                    trialwise_mean=[];
                    trialwise_mean=nanmean(current_siteN1,2);
                    val=[]; BF=[];

                    [val BF]=min(trialwise_mean); % max response has min value

                    BF_responseN1=[];

                    unitwiseBFs(num_sites)=BF;

                    


                    BF_responseN1=zscore(current_siteN1(BF,:));% zscore done here
                    BF_responseP1=zscore(current_siteP1(BF,:));
                    if sum(isnan(current_siteP1(BF,:)))==0  % get rid of nans




                        Pupil_and_Response(unit_cnt).Pupil=PupilSizeFreqXReps(BF,:);
                        Pupil_and_Response(unit_cnt).ResponseN1=BF_responseN1;
                        Pupil_and_Response(unit_cnt).ResponseP1=BF_responseP1;
                        Pupil_and_Response(unit_cnt).walk=WalkBins(BF,:);
                        % RawLFP=[];RawLFP=squeeze(ABR_and_LFP_Data_HS.LFP(fols).RawLFP(num_sites,BF,:,:));
                      



                        trials_with_current_pupil_size=[];
                        for pupbins=1:length(PupilBinEdges)
                            trials_with_current_pupil_size=[];
                            if pupbins==1

                                trials_with_current_pupil_size=find(PupilSizeFreqXReps(BF,:)<=PupilBinEdges(pupbins) & WalkBins(BF,:)==0);
                            else
                                trials_with_current_pupil_size=find(PupilSizeFreqXReps(BF,:)<=PupilBinEdges(pupbins) & PupilSizeFreqXReps(BF,:)>PupilBinEdges(pupbins-1)  & WalkBins(BF,:)==0);
                            end
                            
                            % eval(sprintf('LFPtimeseries(unit_cnt).pupilBin%i=RawLFP(trials_with_current_pupil_size,:);',pupbins))
                            



                            ALL_Binned_data(unit_cnt,pupbins,1)=nanmean(BF_responseN1(trials_with_current_pupil_size)); % mean_response at BF as a function of pupil
                            ALL_Binned_PupilData(unit_cnt,pupbins)=nanmean(PupilSizeFreqXReps(BF,trials_with_current_pupil_size));
                            ALL_sessionwise(unit_session,pupbins,1)= nanmean(BF_responseN1(trials_with_current_pupil_size));

                            ALL_units(unit_cnt).data_trialwise(pupbins).trialsN1=BF_responseN1(trials_with_current_pupil_size);

                            ALL_Binned_data(unit_cnt,pupbins,2)=nanmean(BF_responseP1(trials_with_current_pupil_size)); % mean_response at BF as a function of pupil
                            ALL_sessionwise(unit_session,pupbins,2)= nanmean(BF_responseP1(trials_with_current_pupil_size));

                            ALL_units(unit_cnt).data_trialwise(pupbins).trialsP1=BF_responseP1(trials_with_current_pupil_size);



                        end


                        % add walk bins
                        walk_trials=[];
                        walk_trials=find(WalkBins(BF,:)==1);

                        ALL_Binned_data(unit_cnt,11,1)=mean(BF_responseN1(walk_trials));
                        ALL_sessionwise(unit_session,11,1)=mean(BF_responseN1(walk_trials));
                        ALL_Binned_data(unit_cnt,11,2)=mean(BF_responseP1(walk_trials));
                        ALL_sessionwise(unit_session,11,2)=mean(BF_responseP1(walk_trials));
                        % LFPtimeseries(unit_cnt).pupilBin11=RawLFP(walk_trials,:);




                        unit_cnt=unit_cnt+1;
                        unit_session=unit_session+1;
                    end
                end
            end
        end
    end
    SessionWise_data(sessionCNT)={ALL_sessionwise};
    BFSessionWise(sessionCNT)={unitwiseBFs};
    ALL_sessionwise=[];unitwiseBFs=[];
    sessionCNT = sessionCNT + 1;
end

end



function[Rsqv]=getRsquare_helper(ALL_data,which_peak,ALL_Binned_PupilData)
    Rsqv=[];
    alph=.32;
    % ALL_data=ALL_dataN;
    parfor bootst=1:1000
        get_rand=[];
        data_resampled=[];
        get_rand=randi(size(ALL_data,1),1,size(ALL_data,1));
        data_resampled=ALL_data(get_rand,:);
        data_resampled_xaxis=ALL_Binned_PupilData(get_rand,:);

        inputY_temp=data_resampled(:,1:10); % exclude walking
        inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
        inputX=[];
        inputX_temp=data_resampled_xaxis(:,1:10); % exclude walking
        inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
        % for kk=1:10
        %     inputX=[inputX kk.*ones(1,size(inputY_temp,1))]; % create X axis (independent var)
        % end
        % get R square
        output=[];
        output= fit_quadratic_model_and_cal_R2(inputX, inputY)
        Rsqv(bootst)=output;
    
    end


    if which_peak==-1
        bar(9,mean(Rsqv),.5,'b')
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, Rsqv}, 'type', 'bca','alpha',alph);
        Lci=[];Uci=[];
        Lci=mean(Rsqv)-bca_ci(1);
        Uci=bca_ci(2) - mean(Rsqv);
        errorbar(9,mean(Rsqv),Lci,Uci,'.k')
        N1(1)=mean(Rsqv);
        N1(2)=Lci;
        N1(3)=Uci;

    else
        bar(8,mean(Rsqv),.5,'b')
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, Rsqv}, 'type', 'bca','alpha',alph);
        Lci=[];Uci=[];
        Lci=mean(Rsqv)-bca_ci(1);
        Uci=bca_ci(2) - mean(Rsqv);
        errorbar(8,mean(Rsqv),Lci,Uci,'.k')
         P1(1)=mean(Rsqv);
        P1(2)=Lci;
        P1(3)=Uci;
    end
end


% clearvars -except ABR_and_LFP_Data_HS



function plot_mean_LFPs(ALL_dataN,subP,whichPeak,ALL_Binned_PupilData,SessionWise_data)


pupil_copy=ALL_Binned_PupilData;
for kk=1:length(SessionWise_data)

    temp=size(SessionWise_data{kk},1);

    session_mean_pupil(kk,:)=mean(pupil_copy(1:temp,:));
    pupil_copy(1:temp,:)=[];
end


for kk=1:10
    data=ALL_dataN(:,kk);
    pupilsize=session_mean_pupil(:,kk);

% Number of bootstrap samples
num_bootstrap_samples = 1000;

% Perform BCA bootstrap resampling
% bootstrap_means = bootstrp(num_bootstrap_samples, @mean, data);
bca_confidence_interval(kk,:) = bootci(num_bootstrap_samples, {@mean, data}, 'type', 'bca','alpha',.32);
bca_confidence_interval_xaxis(kk,:) = bootci(num_bootstrap_samples, {@mean,  pupilsize}, 'type', 'bca','alpha',.32);


original_mean(kk) = mean(data);
original_mean_xaxis(kk) = mean(pupilsize);

end
if subP
    

hold on
 errorbar([original_mean_xaxis],original_mean,original_mean'- bca_confidence_interval(:,1) , bca_confidence_interval(:,2) - original_mean',original_mean_xaxis'- bca_confidence_interval_xaxis(:,1) , bca_confidence_interval_xaxis(:,2) - original_mean_xaxis','k');
 plot(original_mean_xaxis,original_mean,'k','LineWidth',2)
else 
    if whichPeak==-1

        
        hold on
        errorbar([original_mean_xaxis],original_mean,original_mean'- bca_confidence_interval(:,1) , bca_confidence_interval(:,2) - original_mean',original_mean_xaxis'- bca_confidence_interval_xaxis(:,1) , bca_confidence_interval_xaxis(:,2) - original_mean_xaxis','m','LineStyle','none','LineWidth',2);
        plot(original_mean_xaxis,original_mean,'m','LineWidth',2)
        errorbar([original_mean_xaxis],original_mean,[], [],original_mean_xaxis'- bca_confidence_interval_xaxis(:,1) , bca_confidence_interval_xaxis(:,2) - original_mean_xaxis','k','LineStyle','none','LineWidth',2);
    else
       
        hold on
        % errorbar([original_mean_xaxis],original_mean,original_mean'- bca_confidence_interval(:,1) , bca_confidence_interval(:,2) - original_mean',original_mean_xaxis'- bca_confidence_interval_xaxis(:,1) , bca_confidence_interval_xaxis(:,2) - original_mean_xaxis','r','LineStyle','none','LineWidth',2);
        errorbar([original_mean_xaxis],original_mean,original_mean'- bca_confidence_interval(:,1) , bca_confidence_interval(:,2) - original_mean',[] , [],'r','LineStyle','none','LineWidth',2);
        
        plot(original_mean_xaxis,original_mean,'r','LineWidth',2)
        % errorbar([original_mean_xaxis],original_mean,[], [],original_mean_xaxis'- bca_confidence_interval_xaxis(:,1) , bca_confidence_interval_xaxis(:,2) - original_mean_xaxis','k','LineStyle','none','LineWidth',2);


    end
end



 parfor bootst=1:1000
    
    get_rand=[];
    data_resampled=[];

    get_rand=randi(size(ALL_dataN,1),1,size(ALL_dataN,1));

    data_resampled=ALL_dataN(get_rand,:);
    data_resampledXaxis=ALL_Binned_PupilData(get_rand,:);

    inputY_temp=data_resampled(:,1:10); % exclude walking
    inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
    inputX=[];
    inputX_temp=data_resampledXaxis(:,1:10);
    inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
    % for kk=1:10
    %     inputX=[inputX kk.*ones(1,size(inputY_temp,1))]; % create X axis (independent var)
    % end
    fitting_coeff=[];
    fitting_coeff=polyfit(inputX,inputY,2);
    y_predicted(:,bootst)=polyval(fitting_coeff,[min(original_mean_xaxis):.01:max(original_mean_xaxis)]);
 end

 if subP
 plot([min(original_mean_xaxis):.01:max(original_mean_xaxis)],mean(y_predicted,2),'r','LineStyle','-.','LineWidth',2)
 end

 % xlim([-1 12])
 ylim([-.1 .1])
 xlabel('frac of max')
 ylabel('Peak amp. (Zscore)')

 set(gca,'TickDir','out')
 set(gca, 'box', 'off')
end




%% plot abr amps

function plot_ABR_amps(ABRdata,subP)

Sessions=ABRdata.Sessions;
     cmap=winter(5);
    ax=[1 2 3 5 6];
    for num_peaks=1:5
    
            peaks_means=[];
            for num_sessions=1:length(Sessions)
                
                
                
                current_level=[];
                current_level=Sessions(num_sessions).sessionNum(ABRdata.ChosenLevel(num_sessions)).Levels;
                
                All_pupilBins_Merged=[];Pupil_bin_ID=[];PupilValue=[];
                for pupbins=1:length(current_level)
                    All_pupilBins_Merged=[All_pupilBins_Merged;current_level(pupbins).peakvalues];
                    Pupil_bin_ID=[Pupil_bin_ID pupbins.*ones(1,length(current_level(pupbins).peakvalues))];
                    PupilValue=[PupilValue current_level(pupbins).pupil_values];
                end
                Zscored_All_pupilBins_Merged=[];
                Zscored_All_pupilBins_Merged=(All_pupilBins_Merged(:,num_peaks));
                
                
                Zscored_All_pupilBins_Merged=zscore(Zscored_All_pupilBins_Merged);
                
                
                for pupbins=1:length(current_level)
                    
                    ABR_in_current_pupil_size=[];
                    ABR_in_current_pupil_size=Zscored_All_pupilBins_Merged(find(Pupil_bin_ID==pupbins));
                    PupilSize_in_current_pupil_size(num_sessions,pupbins)=nanmean(PupilValue(find(Pupil_bin_ID==pupbins)));
                    peaks_means(num_sessions,pupbins)=nanmean(ABR_in_current_pupil_size);
                end
            end


          LALA(num_peaks,:,:) = peaks_means;
            
            for num_bins=1:10

                data_temp=peaks_means(:,num_bins);
                data_temp_xaxis=PupilSize_in_current_pupil_size(:,num_bins);


            bca_ci(num_bins,:) = bootci(1000, {@mean, data_temp}, 'type', 'bca','alpha',.32);
            bca_ci_xaxis(num_bins,:) = bootci(1000, {@mean, data_temp_xaxis}, 'type', 'bca','alpha',.32);
            data_mean(num_bins)=mean(data_temp);
            data_mean_xaxis(num_bins)=mean(data_temp_xaxis);

            end

            if subP
                subplot(2,4,ax(num_peaks))
            end
            
            if subP
                plot([ data_mean_xaxis],data_mean,'k','LineWidth',2)
                hold on
                errorbar([data_mean_xaxis],data_mean,data_mean'-bca_ci(:,1),bca_ci(:,2)-data_mean',data_mean_xaxis'-bca_ci_xaxis(:,1),bca_ci_xaxis(:,2)-data_mean_xaxis','k')
            else
                plot([data_mean_xaxis],data_mean,'Color',cmap(num_peaks,:),'LineWidth',2)
                hold on
                errorbar([data_mean_xaxis],data_mean,data_mean'-bca_ci(:,1),bca_ci(:,2)-data_mean',[],[],'Color',cmap(num_peaks,:),'LineStyle','none','LineWidth',2)
               if num_peaks==5
                errorbar([data_mean_xaxis],data_mean,[],[],data_mean_xaxis'-bca_ci_xaxis(:,1),bca_ci_xaxis(:,2)-data_mean_xaxis','Color','k','LineStyle','none','LineWidth',2)
               end
            end



            y_predicted=[];
            parfor bootst=1:1000

                get_rand=[];
                data_resampled=[];

                get_rand=randi(size(peaks_means,1),1,size(peaks_means,1));

                data_resampled=peaks_means(get_rand,:);
                data_resampled_xaxis=PupilSize_in_current_pupil_size(get_rand,:);

                inputY_temp=data_resampled(:,1:10); % exclude walking
                inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
                inputX=[];
                inputX_temp=data_resampled_xaxis(:,1:10);
                inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
                % for kk=1:10
                %     inputX=[inputX kk.*ones(1,size(inputY_temp,1))]; % create X axis (independent var)
                % end
                fitting_coeff=[];
                fitting_coeff=polyfit(inputX,inputY,2);
                y_predicted(:,bootst)=polyval(fitting_coeff,[min(data_mean_xaxis):.01:max(data_mean_xaxis)]);
            end
            
            if subP
             plot([min(data_mean_xaxis):.01:max(data_mean_xaxis)],mean(y_predicted,2),'r','LineStyle','-.','LineWidth',2)
            end


            ylim([-.1 .1])
            % xlim([-1 12])

            % xticks(2:2:10)
            xlabel('frac of max')

            ylabel('Peak amp. (Zscore)')

            set(gca,'TickDir','out')
            set(gca, 'box', 'off')


      
            
        
    end
    c=2;
end


%% clearvars -except choose_level


function plot_p5p1_ratio(ABRdata,subP)

Sessions=ABRdata.Sessions;
for num_peaks=[1 5]


    if num_peaks==1
        indx=1;
    else indx=2;
    end

    for num_sessions=1:length(Sessions)




        current_level=[];
        current_level=Sessions(num_sessions).sessionNum(ABRdata.ChosenLevel(num_sessions)).Levels;

        All_pupilBins_Merged=[];Pupil_bin_ID=[];PupilValue=[];
        for pupbins=1:length(current_level)
            All_pupilBins_Merged=[All_pupilBins_Merged;current_level(pupbins).peakvalues];
            Pupil_bin_ID=[Pupil_bin_ID pupbins.*ones(1,length(current_level(pupbins).peakvalues))];
            PupilValue=[PupilValue current_level(pupbins).pupil_values];
        end
        Zscored_All_pupilBins_Merged=[];
        Zscored_All_pupilBins_Merged=(All_pupilBins_Merged(:,num_peaks));



        for pupbins=1:length(current_level)

            ABR_in_current_pupil_size=[];
            ABR_in_current_pupil_size=Zscored_All_pupilBins_Merged(find(Pupil_bin_ID==pupbins));
             PupilSize_in_current_pupil_size(num_sessions,pupbins)=nanmean(PupilValue(find(Pupil_bin_ID==pupbins)));
            peaks_means(num_sessions,pupbins,indx)=nanmean(ABR_in_current_pupil_size);
        end
    end


end

    


        p5=squeeze(peaks_means(:,:,2));
        p1=squeeze(peaks_means(:,:,1));
        rat=[];
        rat=p5./p1;
       
        rat(find(rat<0))=NaN;

        log_rat=log(rat);

 

        for num_bins=1:10
            data_temp=log_rat(:,num_bins);
            data_temp_xaxis=PupilSize_in_current_pupil_size(:,num_bins);
            bca_ci(num_bins,:) = bootci(1000, {@nanmean, data_temp}, 'type', 'bca','alpha',.32);
            bca_ci_xaxis(num_bins,:) = bootci(1000, {@nanmean, data_temp_xaxis}, 'type', 'bca','alpha',.32);
            data_mean(num_bins)=nanmean(data_temp);
            data_mean_xaxis(num_bins)=nanmean(data_temp_xaxis);
        end




            % subplot(2,4,7)

            plot([data_mean_xaxis],data_mean,'k','LineWidth',2)
            hold on
            errorbar([data_mean_xaxis],data_mean,data_mean'-bca_ci(:,1),bca_ci(:,2)-data_mean',data_mean_xaxis'-bca_ci_xaxis(:,1),bca_ci_xaxis(:,2)-data_mean_xaxis','k','LineWidth',2,'LineStyle','none')



            y_predicted=[];
            parfor bootst=1:1000

                get_rand=[];
                data_resampled=[];

                get_rand=randi(size(log_rat,1),1,size(log_rat,1));

                data_resampled=log_rat(get_rand,:);
                data_resampled_xaxis=PupilSize_in_current_pupil_size(get_rand,:);

                inputY_temp=data_resampled(:,1:10); % exclude walking
                inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
                inputX=[];
                inputX_temp=data_resampled_xaxis(:,1:10);
                inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
                % for kk=1:10
                %     inputX=[inputX kk.*ones(1,size(inputY_temp,1))]; % create X axis (independent var)
                % end
                fitting_coeff=[];
                t01=find(isnan(inputY)==1);
                inputX(t01)=[];inputY(t01)=[];

                fitting_coeff=polyfit(inputX,inputY,2);
                y_predicted(:,bootst)=polyval(fitting_coeff,[min(data_mean_xaxis):.01:max(data_mean_xaxis)]);
            end
            
             % plot([[min(data_mean_xaxis):.01:max(data_mean_xaxis)]],nanmean(y_predicted,2),'r','LineStyle','-.','LineWidth',2)
            
            ylim([-1.8 -.2])
            % xlim([-1 12])

            % xticks(2:2:10)
            xlabel('frac of max')
            ylabel('log(V/I)')
            set(gca,'TickDir','out')
            set(gca, 'box', 'off')
    
end
      

%% cal R2 with local averaging for LPFs

function[R2_locally_averaged] = get_local_averaged_R2(Pupil_and_Response,whichPeak)

for num_sites=1:length(Pupil_and_Response)

    ylims=[];
    current_site_resp=[];current_site_Pupil=[];
    if whichPeak==-1
        current_site_resp=Pupil_and_Response(num_sites).ResponseN1;
        col='m';
    else
        current_site_resp=Pupil_and_Response(num_sites).ResponseP1;
        col='r';
    end
    current_site_Pupil=Pupil_and_Response(num_sites).Pupil;
    tbd=[];
    tbd=find(Pupil_and_Response(num_sites).walk==1);

    current_site_Pupil(tbd)=[];current_site_resp(tbd)=[];

    sorted_pupil=[];Pupilsorted_resp=[]; ind=[];
    [sorted_pupil ind]=sort(current_site_Pupil);
    Pupilsorted_resp=current_site_resp(ind);

    binning_approach=[1 2 4 6 8];

    for numBA=1:length(binning_approach)
        average_across=[];
        average_across=1:binning_approach(numBA):length(Pupilsorted_resp);
        AV_sorted_pupil=[];
        AV_Pupilsorted_resp=[];

        for num_averages=1:length(average_across)-1

            AV_sorted_pupil(num_averages)=mean(sorted_pupil(average_across(num_averages):average_across(num_averages+1)-1));

            AV_Pupilsorted_resp(num_averages)=mean(Pupilsorted_resp(average_across(num_averages):average_across(num_averages+1)-1));


        end



        numBins(numBA)=length(AV_Pupilsorted_resp);

        sst=[];
        sst = sum((AV_Pupilsorted_resp-mean(AV_Pupilsorted_resp)).^2);
        degree=2;

        dfe=[];dfe=length(AV_Pupilsorted_resp)-3;

        fitting_coeff=[];fitting_coeff=polyfit(AV_sorted_pupil,AV_Pupilsorted_resp,degree);
        y_predicted=[];y_predicted=polyval(fitting_coeff,AV_sorted_pupil);
        SSE_uncv=[];SSE_uncv= sum ( (y_predicted-AV_Pupilsorted_resp).^2);
        rsquare_withoutCV=[];rsquare_withoutCV = 1 - SSE_uncv./sst;
        adjrsquare_withoutCV=[];adjrsquare_withoutCV= 1 - (1-rsquare_withoutCV)*(length(AV_Pupilsorted_resp)-1)/dfe;

        R2_locally_averaged(num_sites,numBA) =adjrsquare_withoutCV;


        [mdl,goodness] = fit(AV_sorted_pupil',AV_Pupilsorted_resp','poly2');
        sse_fit=goodness.sse;
        rsquare_fit=goodness.rsquare;
        R2_locally_averaged(num_sites,numBA)=goodness.adjrsquare;

    end
end

errorbar(mean(R2_locally_averaged),std(R2_locally_averaged)./sqrt(size(R2_locally_averaged,1)),'Color',col,'LineStyle','none','LineWidth',2);
hold on
plot(mean(R2_locally_averaged),'LineWidth',2,'Color',col)
end



%% R2 local averaging for ABR


function[all_comp]=get_R2_ABR_local_averaging(ABRdata);

desired_level=ABRdata.ChosenLevel;
Sessions=ABRdata.Sessions;   
peaks_means=[];
for num_sessions=1:length(Sessions)
    ylims=[];
    for num_peaks=1:5
        current_level=[];
        current_level=Sessions(num_sessions).sessionNum(desired_level(num_sessions)).Levels;
        All_pupilBins_Merged=[];Pupil_bin_ID=[];All_pupilMerged=[];
        for pupbins=1:length(current_level)
            All_pupilBins_Merged=[All_pupilBins_Merged;current_level(pupbins).peakvalues];
            Pupil_bin_ID=[Pupil_bin_ID pupbins.*ones(1,length(current_level(pupbins).peakvalues))];
            All_pupilMerged = [All_pupilMerged current_level(pupbins).pupil_values];
        end
        Zscored_All_pupilBins_Merged=[];
        Zscored_All_pupilBins_Merged=(All_pupilBins_Merged(:,num_peaks));
        raw_peakvalue=[];
        raw_peakvalue=(Zscored_All_pupilBins_Merged);
        Zscored_All_pupilBins_Merged=zscore(Zscored_All_pupilBins_Merged);
        sorted_pupil=[];Pupilsorted_resp=[]; ind=[];
        [sorted_pupil ind]=sort(All_pupilMerged);
        Pupilsorted_resp=Zscored_All_pupilBins_Merged(ind);
        binning_approach=[1 2 4 6 8];
        for numBA=1:length(binning_approach)
            average_across=[];
            average_across=1:binning_approach(numBA):length(Pupilsorted_resp);
            AV_sorted_pupil=[];
            AV_Pupilsorted_resp=[];
            for num_averages=1:length(average_across)-1
                AV_sorted_pupil(num_averages)=mean(sorted_pupil(average_across(num_averages):average_across(num_averages+1)-1));
                AV_Pupilsorted_resp(num_averages)=mean(Pupilsorted_resp(average_across(num_averages):average_across(num_averages+1)-1));
            end
            numBins(numBA)=length(AV_Pupilsorted_resp);
            sst=[];
            sst = sum((AV_Pupilsorted_resp-mean(AV_Pupilsorted_resp)).^2);

            [mdl,goodness] = fit(AV_sorted_pupil',AV_Pupilsorted_resp','poly2');
            sse_fit=goodness.sse;
            rsquare_fit=goodness.rsquare;
            adjrsquare_fit=goodness.adjrsquare;
            all_comp(num_sessions,numBA,num_peaks)=adjrsquare_fit;
        end
    end
    
end
cmap=winter(5);
for numpeak=1:5
    errorbar(mean(all_comp(:,:,numpeak)),std(all_comp(:,:,numpeak))./sqrt(size(all_comp,1)),'Color',cmap(numpeak,:),'LineStyle','none','LineWidth',2)
    plot(mean(all_comp(:,:,numpeak)),'Color',cmap(numpeak,:),'LineWidth',2)
end



end



%% lfp timeseries
function plot_lfp_timeseries_pupil_binned(LFPtimeseries)

cmap1 =flipud(hot(18));
cmap =(winter(10));
 cmap1 = cmap1(8:1:end,:);
 cmap1=[];
 cmap1=cmap;

pre=125; % 125=50 ms, SF=2500

alph=.32; % 68% CI
Smo_Pop_raw=smoothdata(LFPtimeseries,2,'gaussian',25); % 10 ms window smoothing

Resampled_Popdata=resample(Smo_Pop_raw,1,3,'Dimension',2);
loc=.25;
for pupbins=1:10
  
    t=squeeze(Resampled_Popdata(:,:,pupbins));
  
     Lci=[];Uci=[];
    for kk=1:size(t,2)
        temp=t(:,kk);
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, temp}, 'type', 'bca','alpha',alph);
       
        Lci(kk)=mean(temp)-bca_ci(1);
        Uci(kk)=bca_ci(2) - mean(temp);
    end
    % errorbar(1:size(t,2),mean(t),Lci,Uci,'Color',cmap(pupbins,:),'LineWidth',.5)
    % errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap(pupbins,:),'LineWidth',.5,'CapSize',3)
      errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap1(pupbins,:),'LineWidth',.5,'CapSize',3,'LineStyle','none')
    
    hold on
    ax=plot(mean(t(:,25:75)),'Color',cmap(pupbins,:),'LineWidth',1);
     loc=loc+.02;    
    plot([46 49.5],[loc loc],Color=cmap(pupbins,:),LineWidth=1.5)

end
% 
% h1=legend({'','1','','2','','3','','4','','5','','6','','7','','8','','9','','10'},Direction="reverse",Box="off");
% set(h1, 'TextColor','w')


%xlim([25 75])
xticks([1 17.66 51])
xticklabels([-20 0 40])
ylim([-.5 .5])
%rectangle('Position',[26.5,-.41,6,.16],'EdgeColor','g')
rectangle('Position',[26.5,-.402,6,.08],'EdgeColor','g')
rectangle('Position',[20.1,.06,4,.07],'EdgeColor','b')
set(gca, 'box', 'off')
set(gca,'TickDir','out')
xlim([0 52])
xlabel('Time re tone onset (ms)')
ylabel('LFP (zscore)')
text(21.2,.16,'P1','Fontsize',18);
text(28.6,-.425,'N1','Fontsize',18);


  currFig = gcf;
   set(currFig, 'color', 'w');
   annotation('textarrow',[.88 .88],[.735 .89],'FontSize',13,'Linewidth',2)
   t1=text(51.7,.29,'Pupil size','Fontsize',18);
   t1.Rotation=90;
   ax=gca;
   ax.XAxis.FontSize = 18;
   ax.YAxis.FontSize = 18;

insetAx=[0.7338 0.1597 0.1143 0.2858];

axes('Position',insetAx)

for pupbins=1:10
  
    t=squeeze(Resampled_Popdata(:,:,pupbins));
  
     Lci=[];Uci=[];
    for kk=1:size(t,2)
        temp=t(:,kk);
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, temp}, 'type', 'bca','alpha',alph);
       
        Lci(kk)=mean(temp)-bca_ci(1);
        Uci(kk)=bca_ci(2) - mean(temp);
    end
    % errorbar(1:size(t,2),mean(t),Lci,Uci,'Color',cmap(pupbins,:),'LineWidth',.5)
    %errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap(pupbins,:),'LineWidth',.5,'CapSize',3)
      errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap1(pupbins,:),'LineWidth',.5,'CapSize',3,'LineStyle','none')
    hold on
    ax=plot(mean(t(:,25:75)),'Color',cmap(pupbins,:),'LineWidth',1);
    

end
rectangle('Position',[26.5,-.402,6,.07],'EdgeColor','g')



Ax = gca;
Ax.XColor = 'none';
Ax.YColor = 'none';
% ylim([-.41 -.24])
ylim([-.402 -.3320])

xlim([26.5 32.5])

insetAx_P1=[0.2379 0.5820 0.1143 0.2858];
axes('Position',insetAx_P1)
for pupbins=1:10

    t=squeeze(Resampled_Popdata(:,:,pupbins));

    Lci=[];Uci=[];
    for kk=1:size(t,2)
        temp=t(:,kk);
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, temp}, 'type', 'bca','alpha',alph);

        Lci(kk)=mean(temp)-bca_ci(1);
        Uci(kk)=bca_ci(2) - mean(temp);
    end
    % errorbar(1:size(t,2),mean(t),Lci,Uci,'Color',cmap(pupbins,:),'LineWidth',.5)
    % errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap(pupbins,:),'LineWidth',.5)

    errorbar(1:size(t(:,25:75),2),mean(t(:,25:75)),Lci(25:75),Uci(25:75),'Color',cmap1(pupbins,:),'LineWidth',.5,'CapSize',3,'LineStyle','none')
    hold on
    ax=plot(mean(t(:,25:75)),'Color',cmap(pupbins,:),'LineWidth',1);


end
xlim([20.1 24.1])
ylim([.06 .12])

rectangle('Position',[20.1,.06,4,.06],'EdgeColor','b')

Ax = gca;
Ax.XColor = 'none';
Ax.YColor = 'none';

% rectangle('Position',[51,-.48,5,.18],'EdgeColor','g')
% rectangle('Position',[44,.03,4.5,.13],'EdgeColor','b')
set(gca,'TickDir','out')
set(gca, 'box', 'off')

end


%% ABR population timeSeries

function plot_ABR_timeseries(ABRdata)
%cmap =flipud(hot(18));
cmap =(winter(10));
%cmap = cmap(8:1:end,:);

alph=.32; % 68% CI
Smo_Pop_raw=smoothdata(ABRdata,2,'gaussian',15); % 1 ms window smoothing

Resampled_Popdata=resample(Smo_Pop_raw,1,2,'Dimension',2);
loc=.00015;
for pupbins=1:10

    t=squeeze(Resampled_Popdata(:,46:226,pupbins));

    Lci=[];Uci=[];
    for kk=1:size(t,2)
        temp=t(:,kk);
        bca_ci=[];
        bca_ci = bootci(1000, {@mean, temp}, 'type', 'bca','alpha',alph);

        Lci(kk)=mean(temp)-bca_ci(1);
        Uci(kk)=bca_ci(2) - mean(temp);
    end
    % errorbar(1:size(t,2),mean(t),Lci,Uci,'Color',cmap(pupbins,:),'LineWidth',.5)
    errorbar(1:size(t,2),mean(t),Lci,Uci,'Color',cmap(pupbins,:),'LineWidth',.5,'CapSize',1)
    hold on
    ax=plot(mean(t),'Color',cmap(pupbins,:),'LineWidth',1);
    loc=loc+.00001    
    plot([170 180],[loc loc],Color=cmap(pupbins,:),LineWidth=1.5)



end

% legend({'','1','','2','','3','','4','','5','','6','','7','','8','','9','','10'},Direction="reverse",Box="off")
% plot([170 180],[.0015 .0015],Color=cmap(1,:),LineWidth=1.5)
xlim([-10 191])
xticks([1:30:181])
xticklabels([-2:2:10])
ylim([-.0003 .0003])

set(gca,'TickDir','out')
set(gca, 'box', 'off')
ylabel('ABR (mV)')
xlabel('Time re tone onset (ms)')


text(21.2,.16,'P1','Fontsize',18);
text(28.6,-.425,'N1','Fontsize',18);


  currFig = gcf;
   set(currFig, 'color', 'w');
   annotation('textarrow',[.875 .875],[.73 .86],'FontSize',13,'Linewidth',2)
   t1=TextZoomable(187,.000158,'Pupil size','Fontsize',18);
   t1.Rotation=90;
   ax=gca;
   ax.XAxis.FontSize = 18;
   ax.YAxis.FontSize = 18;


end



%% find figure

function h=find_figure(figname)

h = findobj('Tag', figname); % check for pre-existing window
if(isempty(h)); % if none, make one
    h = figure('Tag', figname, 'Name', figname, 'NumberTitle', 'off')
end
figure(h)
end

%% R2 with LOOCV for ABR peaks
function[all_r2s] =getR2_withLOOCV_forABT(ABRdata)
desired_level=ABRdata.ChosenLevel;
Sessions=ABRdata.Sessions;   
peaks_means=[];

alph=.32;
  

for num_peaks=1:5
   
        peaks_means=[];
        for num_sessions=1:length(Sessions)




            current_level=Sessions(num_sessions).sessionNum(desired_level(num_sessions)).Levels;

            All_pupilBins_Merged=[];Pupil_bin_ID=[];PupilValue=[];
            for pupbins=1:length(current_level)
                All_pupilBins_Merged=[All_pupilBins_Merged;current_level(pupbins).peakvalues];
                Pupil_bin_ID=[Pupil_bin_ID pupbins.*ones(1,length(current_level(pupbins).peakvalues))];
                 PupilValue=[PupilValue current_level(pupbins).pupil_values];
            end
            Zscored_All_pupilBins_Merged=[];
            Zscored_All_pupilBins_Merged=(All_pupilBins_Merged(:,num_peaks));


            Zscored_All_pupilBins_Merged=zscore(Zscored_All_pupilBins_Merged);



            for pupbins=1:length(current_level)

                ABR_in_current_pupil_size=[];
                ABR_in_current_pupil_size=Zscored_All_pupilBins_Merged(find(Pupil_bin_ID==pupbins));
                PupilSize_in_current_pupil_size(num_sessions,pupbins)=nanmean(PupilValue(find(Pupil_bin_ID==pupbins)));
                peaks_means(num_sessions,pupbins)=nanmean(ABR_in_current_pupil_size);
            end
        end

        current_ABR_data=[];
        current_ABR_data=peaks_means;

        for bootst=1:1000

            get_rand=[];
            data_resampled=[];inputY_temp=[];inputY=[];

            get_rand=randi(size(current_ABR_data,1),1,size(current_ABR_data,1));

            data_resampled=current_ABR_data(get_rand,:);
            data_resampled_xaxis=PupilSize_in_current_pupil_size(get_rand,:);


            inputY_temp=data_resampled; % exclude walking
            inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
            inputX=[];
            inputX_temp=data_resampled_xaxis; % exclude walking
            inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
            % for kk=1:10
            %     inputX=[inputX kk.*ones(1,size(inputY_temp,1))]; % create X axis (independent var)
            % end

            output=[];
            output= fit_quadratic_model_and_cal_R2(inputX, inputY);


            all_r2s(num_peaks,bootst)=output;
        end

     
    
end






for kk=1:5


        curr=squeeze(all_r2s(kk,:));

        bca_ci = bootci(1000, {@mean, curr}, 'type', 'bca','alpha',alph);

        Lci=mean(curr)-bca_ci(1);
        Uci=bca_ci(2)-mean(curr);

        bar(kk,mean(curr),.5,'b')
        hold on
        errorbar(kk,mean(curr),Lci, Uci,'.k')
        AA(kk,2)=Lci;
        AA(kk,3)=Uci;
        AA(kk,1)=mean(curr);
end

end



%% ratio

function[R2_P5P1_ratio]= getR2_withLOOCV_forABR_peakVI_ratio(ABRdata)
desired_level=ABRdata.ChosenLevel;
Sessions=ABRdata.Sessions;   
peaks_means=[];
alph=.32;

for num_peaks=[1 5]
    if num_peaks==1
        indx=1;
    else indx=2;
    end

        %            peaks_means=[];
        for num_sessions=1:length(Sessions)




            current_level=Sessions(num_sessions).sessionNum(desired_level(num_sessions)).Levels;

            All_pupilBins_Merged=[];Pupil_bin_ID=[];PupilValue=[];
            for pupbins=1:length(current_level)
                All_pupilBins_Merged=[All_pupilBins_Merged;current_level(pupbins).peakvalues];
                Pupil_bin_ID=[Pupil_bin_ID pupbins.*ones(1,length(current_level(pupbins).peakvalues))];
                PupilValue=[PupilValue current_level(pupbins).pupil_values];
            end
            Zscored_All_pupilBins_Merged=[];
            Zscored_All_pupilBins_Merged=(All_pupilBins_Merged(:,num_peaks)); % no zscoring here
           


            for pupbins=1:length(current_level)

                ABR_in_current_pupil_size=[];
                ABR_in_current_pupil_size=Zscored_All_pupilBins_Merged(find(Pupil_bin_ID==pupbins));
                PupilSize_in_current_pupil_size(num_sessions,pupbins)=nanmean(PupilValue(find(Pupil_bin_ID==pupbins)));
                peaks_means(num_sessions,pupbins,indx)=nanmean(ABR_in_current_pupil_size);
            end
        end


        

 
end

    
    
        
        p5=squeeze(peaks_means(:,:,2));
        p1=squeeze(peaks_means(:,:,1));
        rat=[];
        rat=p5./p1;
       

        rat(find(rat<0))=NaN; % remove negative values before taking log

        log_rat=log(rat);

          for bootst=1:1000

            get_rand=[];
            data_resampled=[];

            get_rand=randi(size(log_rat,1),1,size(log_rat,1));

            data_resampled=log_rat(get_rand,:);
            data_resampled_xaxis=PupilSize_in_current_pupil_size(get_rand,:);


            inputY_temp=data_resampled; % exclude walking
            inputY=reshape(inputY_temp,1,prod(size(inputY_temp))); % flatten the data
            inputX=[];
            inputX_temp=data_resampled_xaxis; % exclude walking
            inputX=reshape(inputX_temp,1,prod(size(inputX_temp)));
                t01=find(isnan(inputY)==1);
                inputX(t01)=[];inputY(t01)=[];

            output=[];
            output= fit_quadratic_model_and_cal_R2(inputX, inputY);

            R2_P5P1_ratio(bootst)=output;
          end


          bca_ci = bootci(1000, {@nanmean, R2_P5P1_ratio}, 'type', 'bca','alpha',alph);

          Lci=nanmean(R2_P5P1_ratio)-bca_ci(1);
          Uci=bca_ci(2)-nanmean(R2_P5P1_ratio);

          bar(6,nanmean(R2_P5P1_ratio),.5,'b')
          hold on
          errorbar(6,nanmean(R2_P5P1_ratio),Lci, Uci,'.k')

end



%% V5

function R2sq= fit_quadratic_model_and_cal_R2(inputX, inputY)

% fit quadratic model and calculate explained variance with LOO cross
% validation


degree=2;

parfor num_cv=1:length(inputY)
    inputY_copy=[]; inputX_copy=[];
    inputY_copy=inputY;
    inputX_copy=inputX;

    testY=inputY_copy(num_cv);
    testX=inputX_copy(num_cv);

    inputY_copy(num_cv)=[];
    inputX_copy(num_cv)=[];


    fitting_coeff=polyfit(inputX_copy,inputY_copy,degree);
    y_predicted=polyval(fitting_coeff,testX);

    SQV_temp(num_cv) = (testY-y_predicted).^2;

    PRED(num_cv)=y_predicted;
    YY(num_cv)=testY;
end


SQV=sum(SQV_temp);

meanY=mean(inputY);

SST = sum((meanY-inputY).^2);

R2sq = 1 - SQV./SST;


end

%% text zoomable

function txtHandle = TextZoomable(x,y,varargin)
% txtHandle = FixedSizeText(x,y,varargin)
%
% Version 2.0, 14 July 2015
% Adds text to a figure in the normal manner, except that this text
% grows/shrinks with figure scaling and zooming, unlike normal text that
% stays at a fixed font size during figure operations. Note it scales with
% figure height - for best scaling use 'axis equal' before setting up the
% text.
%
% All varargin{:} arguments will be passed directly on to the text
% function (text properties, etc.)
%
% (doesn't behave well with FontUnits = 'normalized')
%
% example:
%
% figure(1); clf;
% rectangle('Position', [0 0 1 1]);
% rectangle('Position', [.25 .25 .5 .5]);
% 
% th = TextZoomable(.5, .5, 'red', 'color', [1 0 0], 'Clipping', 'on');
% th2 = TextZoomable(.5, .1, 'blue', 'color', [0 0 1]);
%
%
% Ken Purchase, 4-25-2013, with many thanks to the Matlab User
% Community, including Matt J, Hoi Wong, Philip Caplan.
%
% Modified by Philip Caplan on 12/10/2014 to support R2014b
%
% Modified by Hoi Wong on 1/23/2014 so that the 'UserData' will not conflict
% with what is already set by gscatter() if previously called.
% create the text
%
% Modified by Ken, incorporating above changes and version compatibility
% for previous versions
%
% Note there is an issue on my computer that I haven't yet fixed -
% initially the font displays in an absurdly wrong size, but as soon as I
% zoom in or out, all is correct. I don't have time to dig into it now, but
% if you fix it, give me a comment on the TextZoomable page at the file
% exchange, and I can publish your fix. Many thanks.
%
%    
    
    txtHandle = text(x,y,varargin{:});
    % NOTE: The function can be disabled by returning here.
    
    % detect its size relative to the figure, and set up listeners to resize
    % it as the figure resizes, or axis limits are changed.
    hAx = gca;
    hFig = get(hAx,'Parent');
    
    fs = get(txtHandle, 'FontSize');
    ratios = fs * diff(get(hAx,'YLim')) / max(get(hFig,'Position') .* [0 0 0 1]);
    
    % append the handles and ratios to the user data - repeated calls will
    % add each block of text to the list
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(hAx, 'TextZoomable_UserData');
    else
      ud = get(hAx, 'UserData');
    end
    
    if isfield(ud, 'ratios')
      ud.ratios = [ud.ratios(:); ratios];
    else
      ud.ratios = ratios;
    end
    if isfield(ud, 'handles')
      ud.handles = [ud.handles(:); txtHandle];
    else
      ud.handles = txtHandle;
    end
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      setappdata(hAx, 'TextZoomable_UserData', ud);
    else
      set(hAx,'UserData', ud);
    end
    
    localSetupPositionListener(hFig,hAx);
    localSetupLimitListener(hAx);
end
  

%% Helper Functions
  function fs = getBestFontSize(imAxes)
    % Try to keep font size reasonable for text
    hFig = get(imAxes,'Parent');
    hFigFactor = max(get(hFig,'Position') .* [0 0 0 1]);
    axHeight = diff(get(imAxes,'YLim'));
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');  % stored in teh first user data.
    end
    
    fs = round(ud.ratios * hFigFactor / axHeight);
    fs = max(fs, 3);
  end
  
  
  function localSetupPositionListener(hFig,imAxes)
    % helper function to sets up listeners for resizing, so we can detect if
    % we would need to change the fontsize
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      PostPositionListener = addlistener(hFig,'SizeChanged',...
        @(o,e) localPostPositionListener(o,e,imAxes) );
    else
      PostPositionListener = handle.listener(hFig,'ResizeEvent',...
        {@localPostPositionListener,imAxes});
    end
    
    setappdata(hFig,'KenFigResizeListeners',PostPositionListener);
  end
  
  
  function localPostPositionListener(~,~,imAxes)
    % when called, rescale all fonts in image
    if( ~isvalid(imAxes) ) % The imAxes might be deleted, which getappdata() doesn't have an output and will throw an exception
      return;
    end
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');
    end
    
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
      set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end
  end
  
  
  function localSetupLimitListener(imAxes)
    % helper function to sets up listeners for zooming, so we can detect if
    % we would need to change the fontsize
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      LimListener = addlistener(imAxes,{'XLim','YLim'},'PostSet',@localLimitListener);
    else
      hgp     = findpackage('hg');
      axesC   = findclass(hgp,'axes');
      LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
        'PropertyPostSet',@localLimitListener);      
    end
    
    hFig = get(imAxes,'Parent');
    setappdata(hFig,'KenAxeResizeListeners',LimListener);
  end
  
  
  function localLimitListener(~,event)
    % when called, rescale all fonts in image
    imAxes = event.AffectedObject;
    
    if ~verLessThan('matlab', '8.4.0') % R2014b or later.
      ud = getappdata(imAxes, 'TextZoomable_UserData');
    else
      ud = get(imAxes,'UserData');
    end
    
    fs = getBestFontSize(imAxes);
    for ii = 1:length(ud.handles)
      set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
    end
  end
%pc



























% 
% 
% 
%     % create the text
%     txtHandle = text(x,y,varargin{:});
% 
%     % detect its size relative to the figure, and set up listeners to resize
%     % it as the figure resizes, or axis limits are changed.
%     hAx = gca;
%     hFig = get(hAx,'Parent');
%     
%     fs = get(txtHandle, 'FontSize');
%     ratios = fs * diff(get(hAx,'YLim')) / max(get(hFig,'Position') .* [0 0 0 1]);
%     
%     
%     % append the handles and ratios to the user data - repeated calls will
%     % add each block of text to the list
%     ud = get(hAx, 'UserData');
%     if isfield(ud, 'ratios')
%         ud.ratios = [ud.ratios(:); ratios];
%     else
%         ud.ratios = ratios;
%     end
%     if isfield(ud, 'handles')
%         ud.handles = [ud.handles(:); txtHandle];
%     else
%         ud.handles = txtHandle;
%     end
%     
%     set(hAx,'UserData', ud);
%     localSetupPositionListener(hFig,hAx);
%     localSetupLimitListener(hAx);
% 
% end
% 
% 
% %% Helper Functions
% 
% function fs = getBestFontSize(imAxes)
%     % Try to keep font size reasonable for text
%     hFig = get(imAxes,'Parent');
%     hFigFactor = max(get(hFig,'Position') .* [0 0 0 1]);  
%     axHeight = diff(get(imAxes,'YLim'));
%     ud = get(imAxes,'UserData');  % stored in teh first user data.
%     fs = round(ud.ratios * hFigFactor / axHeight);    
%     fs = max(fs, 3);
% end
% 
% function localSetupPositionListener(hFig,imAxes)
%     % helper function to sets up listeners for resizing, so we can detect if
%     % we would need to change the fontsize
%     PostPositionListener = handle.listener(hFig,'ResizeEvent',...
%         {@localPostPositionListener,imAxes});
%     setappdata(hFig,'KenFigResizeListeners',PostPositionListener);
% end
% 
% function localPostPositionListener(~,~,imAxes) 
%     % when called, rescale all fonts in image
%     ud = get(imAxes,'UserData');
%     fs = getBestFontSize(imAxes);
%     for ii = 1:length(ud.handles)
%         set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
%     end   
% end
% 
% function localSetupLimitListener(imAxes)
%     % helper function to sets up listeners for zooming, so we can detect if
%     % we would need to change the fontsize
%     hgp     = findpackage('hg');
%     axesC   = findclass(hgp,'axes');
%     LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
%         'PropertyPostSet',@localLimitListener);
%     hFig = get(imAxes,'Parent');
%     setappdata(hFig,'KenAxeResizeListeners',LimListener);
% end
% 
% function localLimitListener(~,event)
%     % when called, rescale all fonts in image
%     imAxes = event.AffectedObject;
%     ud = get(imAxes,'UserData');
%     fs = getBestFontSize(imAxes);
%     for ii = 1:length(ud.handles)
%         set(ud.handles(ii),'fontsize',fs(ii),'visible','on');
%     end
% end

