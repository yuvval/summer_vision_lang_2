function [seq, score_track, score_history] = viterbi_yuval(em_scores, tr_scores, last_frame_states_mask, prev_frame_scores, t)

Nframes = length(em_scores);

n_states = length(em_scores{t});

if t>1
%     prev_frame_scores
    curr_tr_scores = tr_scores{t-1};
else    
    curr_tr_scores = -1*ones(1, n_states); % uniform prior distribution
end

frame_scores = nan(1, n_states);
best_income_transitions = nan(1, n_states);
for s=1:n_states
    [mx_score, best_income_transitions(s)] = max(prev_frame_scores(:) + curr_tr_scores(:,s));
    frame_scores(s) = em_scores{t}(s) + mx_score;
end

% [best_score, best_state] = max(frame_scores); % A BUG!!

if t < Nframes
    [seq, score_track, score_history] = viterbi_yuval(em_scores, tr_scores, last_frame_states_mask, frame_scores, t+1);
end

if t == Nframes
    score_history = cell(t,1);
    [seq, score_track] = deal(nan(t+1,1));
    frame_scores(~last_frame_states_mask) = -inf;
    [~, best_state] = max(frame_scores);
    seq(end-1:end) = [best_income_transitions(best_state); best_state];
    score_history{t} = frame_scores;
    score_track(t+1) = frame_scores(best_state);
%     frames_scores = frame_scores;
else
    seq(t) = best_income_transitions(seq(t+1));
    score_track(t+1) = frame_scores(seq(t+1));
    score_history{t} = frame_scores;
%     frames_scores = [frame_scores; frames_scores];
end

if t == 1
    % trim init state
    seq(1) = [];
    score_track(1) = [];    
end

