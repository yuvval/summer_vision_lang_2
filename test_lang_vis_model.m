close all
clear


%% 2D
%     ppvid = load('preprocessed_videos/2chairs_approach_diagonal_detections_thm1_05_top_7.mat');
%     ppvid = load('preprocessed_videos/2chairs_approach_behind_detections_thm1_05_top_7.mat');
%     ppvid = load('preprocessed_videos/approach_behind_detections_thm1_05_top_7.mat');
%     ppvid = load('preprocessed_videos/2chairs_approach_side_detections_thm1_05_top_7.mat');
%      ppvid = load('preprocessed_videos/approach_people_side_detections_thm1_05_top_7.mat');
%      ppvid = load('preprocessed_videos/approach_people_behind_detections_thm1_05_top_7.mat');
%      ppvid = load('preבג processed_videos/approach_people_diagonal_detections_thm1_05_top_7.mat');
%      ppvid = load('preprocessed_videos/after2sec_behind_detections_thm1_5_top_3.mat');
     ppvid = load('preprocessed_videos/after2sec_behind_detections_thm1_5_top_7.mat');
%      ppvid = load('preprocessed_videos/after2sec_diagonal_detections_thm1_5_top_3.mat');

%% 3D
%     ppvid = load('preprocessed_videos/2chairs_approach_diagonal_3D_detections_thm1_05_top_7.mat');
%     ppvid = load('preprocessed_videos/2chairs_approach_side_3D_detections_thm1_05_top_7.mat');
%    ppvid = load('preprocessed_videos/2chairs_approach_behind_3D_detections_thm1_5_top_7.mat');
%     ppvid = load('preprocessed_videos/approach_people_behind_3D_detections_thm1_05_top_7.mat');
%    ppvid = load('preprocessed_videos/after2sec_behind_3D_detections_thm1_5_top_3.mat');
%      ppvid = load('preprocessed_videos/after2sec_diagonal_3D_detections_thm1_5_top_3.mat');

    % setting the tuning params for probabilities and features binning / sigmoiding
    % emission probablities sigmoid params
    tuning_params.other.sig_a = 10;
    tuning_params.other.sig_b = -0.8;
    
    tuning_params.person.sig_a = 5;
    tuning_params.person.sig_b = -0.4;
    
    tuning_params.chair.sig_a = 10;
    tuning_params.chair.sig_b = -0.87;
    
    % transition probablities sigmoid params
    tuning_params.sig_a_trans = 0.3;
    tuning_params.sig_b_trans = -4;
    
    [tracker_scores.em, tracker_scores.tr, tracker_feats] = generate_scores_from_2d_preprocessed_video(ppvid, tuning_params);
    
%     verb = 'is on right side of';
%     noun1 = 'chair';
    verb = 'approach';
    noun1 = 'person';

    noun2 = 'chair';
    [ cross_em_scores, cross_tr_scores_mat, cross_p_all_hmms_states, debug_info ] = eval_cross_prod_trellis( verb, noun1, noun2, tracker_scores, tracker_feats);
    
    last_frame_states_mask = cross_p_all_hmms_states{end}(3,:) == 3;
    [seq, score_track, score_history] = viterbi_yuval(cross_em_scores, cross_tr_scores_mat, last_frame_states_mask, 0, 1);
    
%% visualize sequence
if true
    frame_sample_interval = 15;
    obj = VideoReader(['voc-dpm/' ppvid.vid_fname]);
    video = obj.read();
    
    boxes = ppvid.boxes;
    
    t=1;
    feat_history = [];
    trk_history = nan(5,1);
%     for k = 1:frame_sample_interval:size(video,4)
trim_first_seconds = 3;
for k=1 + (trim_first_seconds*30):frame_sample_interval:size(video,4)
        
        im=video(:,:,:,k);
        imshow(im);
        
        if t>1
            d_prev = seq(t-1);
        else
            d_prev = 1;
        end
        
        trk_history(:, end+1) = cross_p_all_hmms_states{t}(:, seq(t));
        
        colors = {'r', 'b'};
        for trkr = 1:2
            
            
            d = cross_p_all_hmms_states{t}(trkr, seq(t));
            x1 = boxes{t}(d,1);
            x2 = boxes{t}(d,2);
            y1 = boxes{t}(d,3);
            y2 = boxes{t}(d,4);
            label = ppvid.classes_names{ppvid.classes{t}(d)};
            %     label = sprintf('%s, %2.3f', label, ppvid.scores{t}(d));
            feat_name = 'velocity_abs';
            feat_id = find(ismember(tracker_feats.names, feat_name));
            
            %     feat_val = ppvid.scores{t}(d);
%             feat_val = tracker_feats.values{t}(d_prev, d, feat_id);
            feat_val = score_track(t);
            feat_history(end+1) = feat_val;
            label = sprintf('%s, %2.3f', label, feat_val);
            line([x1 x1 x2 x2 x1]', [y1 y2 y2 y1 y1]', 'color', colors{trkr}, 'linewidth', 3, 'linestyle', '-');
            text(x1, y1, label, 'Color', 'white');
        end
        drawnow;
        shg;
        pause(0.9)
        t=t+1;
        
    end
end

trk_history

% figure
% hist(feat_history, 0:0.25:5);shg