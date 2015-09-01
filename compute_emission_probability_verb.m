%%(Pv=stop, Cv=stop, rel.dist!=close)

function [em_prob_verb] = compute_emission_probability_verb(verb_name, verb_state_n, tracker_feats, n_frame, n_tracker1_state, n_tracker2_state)
%tracker_feats=feat_per_tr;
%I should find if the class of each state (detection) of each frame is a certain noun

switch verb_name
    case 'is on right side of' 
        em_prob_verb = right_side_of_em_prob(verb_state_n, tracker_feats, n_frame, n_tracker1_state, n_tracker2_state);
    case 'approach' 
        em_prob_verb = approach_em_prob(verb_state_n, tracker_feats, n_frame, n_tracker1_state, n_tracker2_state);
end
end

function [em_prob_verb] = right_side_of_em_prob(verb_state_n, tracker_feats, n_frame, n_tracker1_state, n_tracker2_state)
feat_name = 'velocity_binned';
feat_id = find(ismember(tracker_feats.names, feat_name));
velocity1_binned = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
velocity2_binned = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

feat_name = 'center_x';
feat_id = find(ismember(tracker_feats.names, feat_name));
center1_x = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
center2_x = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

feat_name = 'center_y';
feat_id = find(ismember(tracker_feats.names, feat_name));
center1_y = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
center2_y = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

% distance= sqrt((center1_x-center2_x)^2+(center1_y-center2_y)^2);

my_eps = 1e-200;
% my_eps = 0;
if verb_state_n == 1
%     if velocity1_binned==1 && velocity2_binned==1 &&  (center1_x > center2_x)
    if  (center1_x > center2_x)
        em_prob_verb = 1;
    else 
        em_prob_verb = my_eps;
    end
end

end

function [em_prob_verb] = approach_em_prob(verb_state_n, tracker_feats, n_frame, n_tracker1_state, n_tracker2_state)
feat_name = 'velocity_binned';
feat_id = find(ismember(tracker_feats.names, feat_name));
velocity1_binned = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
velocity2_binned = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

feat_name = 'center_x';
feat_id = find(ismember(tracker_feats.names, feat_name));
center1_x = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
center2_x = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

feat_name = 'center_y';
feat_id = find(ismember(tracker_feats.names, feat_name));
center1_y = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
center2_y = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);

if any(ismember(tracker_feats.names, 'center_z'))
    % With depth
    % center z coordinate
    feat_name = 'center_z';
    feat_id = find(ismember(tracker_feats.names, feat_name));
    center1_z = tracker_feats.values{n_frame}(1, n_tracker1_state, feat_id);
    center2_z = tracker_feats.values{n_frame}(1, n_tracker2_state, feat_id);
else
    center1_z = 0;
    center2_z = 0;
end
distance= sqrt((center1_x-center2_x)^2+(center1_y-center2_y)^2 + (center1_z-center2_z)^2);

my_eps = 1e-200;
% my_eps = 0;
if verb_state_n == 1
    if velocity1_binned==1 && velocity2_binned==1 &&  (distance>100)
        em_prob_verb = 1;
    else 
        em_prob_verb = my_eps;
    end
end
        
if verb_state_n == 2
%     if velocity1_binned==2 && velocity2_binned==1 &&  (distance>100)
    if (distance>100)
        em_prob_verb = 1;
    else 
        em_prob_verb = my_eps;
    end
end

                      
if verb_state_n == 3
%     if velocity1_binned==1 && velocity2_binned==1 &&  (distance < 100) &&(distance > 1)  && n_frame>1 
    if (distance < 100) &&(distance > 1)  && n_frame>1 
        em_prob_verb = 1;
    else 
        em_prob_verb = my_eps;
    end
end

end