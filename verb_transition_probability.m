function verb_tr_prob = verb_transition_probability(verb)

switch verb
    case 'approaches'
        verb_tr_prob = [0.9 0.1 0;0 0.9 0.1; 0 0 1];
%         verb_tr_prob = [0.9 0.1 0;0.0 0.9 0.1; 0 0.01 0.99];
    case 'is on right side of' 
        verb_tr_prob = [1];
    otherwise
        error('unknown verb')
end

        