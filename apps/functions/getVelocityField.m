function v = getVelocityField(x, task)

    scale = 0.5;

    % v
    inv_x = locally_weighted_translation_inverse(task.rho,task.c,task.v,x);
    J = locally_weighted_translation_derivative(task.rho,task.c,task.v,inv_x);
    v = J * task.A * (inv_x - task.qf);

    if norm(v) > scale 
        v = scale * v/norm(v);
    end
    
end