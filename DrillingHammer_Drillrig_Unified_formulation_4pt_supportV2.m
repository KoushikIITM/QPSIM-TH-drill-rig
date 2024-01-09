clc;
clear;

% This code is to simulate a boom mounted drilling hammer


function rx=Rot_x(x_angle)
  rx= [1 0 0; 0 cos(x_angle) -sin(x_angle); 0 sin(x_angle) cos(x_angle)];
endfunction

function ry=Rot_y(y_angle)
  ry= [cos(y_angle) 0 sin(y_angle); 0 1 0; -sin(y_angle) 0 cos(y_angle)];
endfunction

function rz=Rot_z(z_angle)
  rz= [cos(z_angle) -sin(z_angle) 0; sin(z_angle) cos(z_angle) 0; 0 0 1];
endfunction

function ret_a = Ai_mat(phi,theta,psi)
  ret_a=Rot_z(phi)*Rot_x(theta)*Rot_z(psi);
endfunction

function ret_g = Gi_mat(phi, theta, psi)
  ret_g=[0 cos(phi) sin(theta)*sin(phi);...
         0 sin(phi) -sin(theta)*cos(phi);...
         1 0 cos(theta)];
endfunction

function skew_mat = skew_sym(mat)
  skew_mat=[0 -mat(3) mat(2);...
            mat(3) 0 -mat(1);...
            -mat(2) mat(1) 0];
endfunction

function ret_mat = GiThitadoti_thetai(phi,theta,psi, phi_dot, theta_dot, psi_dot) %derivative w.r.t. Euler angles of G^i \dot \theta
ret_mat(1,:) = [-theta_dot*sin(phi)+psi_dot*sin(theta)*cos(phi) psi_dot*cos(theta)*sin(psi) 0];
ret_mat(2,:) = [theta_dot*cos(phi)+psi_dot*sin(theta)*sin(phi) -psi_dot*cos(theta)*cos(psi) 0];
ret_mat(3,:) = [0 -psi_dot*sin(theta) 0];
endfunction

function ret_mat = GibarThitadoti_thetai(phi,theta,psi, phi_dot, theta_dot, psi_dot) %derivative w.r.t. Euler angles of \bar G^i \dot \theta
ret_mat(1,:) = [0 phi_dot*cos(theta)*sin(psi) phi_dot*sin(theta)*cos(psi)-theta_dot*sin(psi)];
ret_mat(2,:) = [0 phi_dot*cos(theta)*cos(psi) -phi_dot*sin(theta)*sin(psi)-theta_dot*cos(psi) ];
ret_mat(3,:) = [0 -phi_dot*sin(theta) 0];
endfunction

function ret_mat = dHqidot_dqi(A_i,G_i,G_i_bar,u_i_bar,q_i,q_i_dot) %function to return the Jacobian of H matrix w.r.t q_i
 ret_mat=[zeros(3,3) skew_sym(A_i*skew_sym(u_i_bar)'*G_i_bar*[q_i_dot(4);q_i_dot(5);q_i_dot(6)])'*G_i+A_i*skew_sym(u_i_bar)'*GibarThitadoti_thetai(q_i(4),q_i(5),q_i(6),q_i_dot(4),q_i_dot(5),q_i_dot(6))];
 %skew_sym(A_i*skew_sym(u_i_bar)'*G_i_bar*[q_i_dot(4);q_i_dot(5);q_i_dot(6)])'*G_i
 %A_i*skew_sym(skew_sym(u_i_bar)'*G_i_bar*[q_i_dot(4);q_i_dot(5);q_i_dot(6)])'*G_i_bar
 %A_i*skew_sym(u_i_bar)'*GibarThitadoti_thetai(q_i(4),q_i(5),q_i(6),q_i_dot(4),q_i_dot(5),q_i_dot(6))
 %A_i*skew_sym(u_i_bar)'*GibarThitadoti_thetai(q_i(4),q_i(5),q_i(6),q_i_dot(4),q_i_dot(5),q_i_dot(6))
 endfunction
                     
function ret_mat = dvi_dqi(A_i,G_i,v_i_bar)
  ret_mat=[zeros(3,3) skew_sym(A_i*v_i_bar)'*G_i];
endfunction

function ret_mat = drij_dqi(A_i,G_i,u_i_p_bar)
  ret_mat=[eye(3,3) skew_sym(A_i*u_i_p_bar)'*G_i];
endfunction

function ret_mat = drij_dqj(A_j,G_j,u_j_p_bar)
  ret_mat=[-eye(3,3) -skew_sym(A_j*u_j_p_bar)'*G_j];
endfunction

function ret_mat = calc_M_i(link_i)
  ret_mat=zeros(6,6);
  ret_mat(1:3,1:3)=link_i.mass*eye(3,3);
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i_bar=A_i'*Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  ret_mat(4:6,4:6)= G_i_bar'*link_i.I_theta_bar*G_i_bar;
endfunction

function ret_mat = Gi_bar_dot(phi, theta, psi, phi_dot, theta_dot, psi_dot)
 ret_mat(1,1:3)=[theta_dot*cos(theta)*sin(psi)+psi_dot*sin(theta)*cos(psi), -psi_dot*sin(psi), 0];
 ret_mat(2,1:3)=[theta_dot*cos(theta)*cos(psi)-psi_dot*sin(theta)*sin(psi), -psi_dot*cos(psi), 0];
 ret_mat(3,1:3)=[-theta_dot*sin(theta), 0, 0];
endfunction

function ret_mat = Qv_theta(link_i)
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i_bar=A_i'*Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i_bar_dot=Gi_bar_dot(link_i.cm(4),link_i.cm(5),link_i.cm(6),link_i.vel(4),link_i.vel(5),link_i.vel(6));
  omega_bar=G_i_bar_dot*[link_i.vel(4);link_i.vel(5);link_i.vel(6)]; %%CHECK
  Qv_theta=-G_i_bar'*[skew_sym(omega_bar)*link_i.I_theta_bar*omega_bar+link_i.I_theta_bar*G_i_bar_dot*[link_i.vel(4);link_i.vel(5);link_i.vel(6)]];
  global Qv_sys;
  Qv_sys(link_i.body_index*6-2:link_i.body_index*6,1)=Qv_theta;
endfunction

function ret_mat= calc_Qe(force_i,link_i,pt)
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));

  global Qe_sys;
  Qe_sys(link_i.body_index*6-5:link_i.body_index*6-3,1)=Qe_sys(link_i.body_index*6-5:link_i.body_index*6-3,1)+force_i.val;
  Qe_sys(link_i.body_index*6-2:link_i.body_index*6,1)=Qe_sys(link_i.body_index*6-2:link_i.body_index*6,1)-G_i'*skew_sym(A_i*pt)'*force_i.val;
  %Qe_sys(link_i.body_index*6-2:link_i.body_index*6,1)
endfunction

function ret_mat= J_mbd(link_i,pt,pt_index)
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  body_index=link_i.body_index;
  global J_mbd_sys;
  J_mbd_sys(body_index*6-5:body_index*6-3,pt_index*3-2:pt_index*3)= eye(3,3);
  J_mbd_sys(body_index*6-2:body_index*6,pt_index*3-2:pt_index*3)=-G_i'*skew_sym(A_i*pt)';
endfunction 

function ret_mat= Li_mat(link_i,pt,pt_index)
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  body_index=link_i.body_index;
  global L_i_sys;
  L_i_sys(pt_index*3-2:pt_index*3,body_index*6-5:body_index*6-3)=eye(3,3);
  L_i_sys(pt_index*3-2:pt_index*3,body_index*6-2:body_index*6)=-skew_sym(A_i*pt)*G_i;
endfunction


function ret_joint_data = spherical_joint(link_i, pt_i,joint_axis_i, link_j, pt_j,joint_axis_j )
  %joint_axis_i is in BCS of body i. joint_axis_i(:,1), joint_axis_i(:,2) and joint_axis_i(:,3) form an orthogonal triad 
  %joint_axis_j is in BCS of body j. joint_axis_j(:,2) and joint_axis_i(:,2) are perpendicular to each other
  %joint_axis_j(:,2) and joint_axis_i(:,2) are always perpendicular to each other restricting relative rotation
  
  %eq. 7.31 from Computational dynamics by AA Shabana
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  A_j=Ai_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  %eq. 7.72 from Computational dynamics by AA Shabana
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_j=Gi_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  
  G_i_bar=A_i'*G_i;
  G_j_bar=A_j'*G_j;
  
  %eq. 7.197 from Computational dynamics by AA Shabana
  C= link_i.cm(1:3)+A_i*pt_i - (link_j.cm(1:3)+A_j*pt_j);
  
  H_i_p=[eye(3) skew_sym(A_i*pt_i)'*G_i]; %eq. 7.205 from Computational dynamics by AA Shabana
  H_j_p=[eye(3) skew_sym(A_j*pt_j)'*G_j]; %eq. 7.205 from Computational dynamics by AA Shabana

  
  %eq. 7.204 from Computational dynamics by AA Shabana
  Cq=[H_i_p -H_j_p];
  
  ret_joint_data.C=C;
  ret_joint_data.Cq=Cq;
  ret_joint_data.nc=3; % no. of constraints
  
  global Nc;
  global C_sys;
  global Cq_sys;
  global Cq_sys;
  global joint_index;
  global C_sys_pointer;
  global Qd_sys;
  global Cq_sys_dyn;
  joint_index=joint_index+1;

  C_sys(C_sys_pointer:C_sys_pointer+2,1)=C;
  Cq_sys(C_sys_pointer:C_sys_pointer+2,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys(C_sys_pointer:C_sys_pointer+2,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+2,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+2,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  
  q_i_dot=link_i.vel;
  q_j_dot=link_j.vel;
  q_i=link_i.cm;
  q_j=link_j.cm;
  
  %Qd=[zeros(3,3)  skew_sym(A_i*skew_sym(pt_i)'*G_i_bar*[q_i_dot(4);q_i_dot(5);q_i_dot(6)])'*G_i+A_i*skew_sym(pt_i)'*GibarThitadoti_thetai(link_i.cm(4),link_i.cm(5),link_i.cm(6),q_i_dot(4),q_i_dot(5),q_i_dot(6)) ...
  %   -zeros(3,3) -skew_sym(A_j*skew_sym(pt_j)'*G_j_bar*[q_j_dot(4);q_j_dot(5);q_j_dot(6)])'*G_j-A_j*skew_sym(pt_j)'*GibarThitadoti_thetai(link_j.cm(4),link_j.cm(5),link_j.cm(6),q_j_dot(4),q_j_dot(5),q_j_dot(6))]*[q_i_dot;q_j_dot]
  
  Qd=-[dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot) -dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)]*[q_i_dot;q_j_dot];
  
  Qd_sys(C_sys_pointer:C_sys_pointer+2,1) = Qd;
  C_sys_pointer=C_sys_pointer+3;
  %Nc=nc;
endfunction

function ret_joint_data = cylindrical_joint(link_i, pt_i,joint_axis_i, link_j, pt_j,joint_axis_j )
  %joint_axis_i is in BCS of body i. joint_axis_i(:,1), joint_axis_i(:,2) and joint_axis_i(:,3) form an orthogonal triad 
  %joint_axis_j is in BCS of body j
  
  %eq. 7.31 from Computational dynamics by AA Shabana
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  A_j=Ai_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  %eq. 7.72 from Computational dynamics by AA Shabana
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_j=Gi_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  
  G_i_bar=A_i'*G_i;
  G_j_bar=A_j'*G_j;
  
  r_ij=link_i.cm(1:3)+A_i*pt_i - (link_j.cm(1:3)+A_j*pt_j);
  v_i_1= A_i*joint_axis_i(:,2);
  v_i_2= A_i*joint_axis_i(:,3);
  v_j=A_j*joint_axis_j(:,1);
  v_i_1_bar= joint_axis_i(:,2);
  v_i_2_bar= joint_axis_i(:,3);
  v_j_bar=joint_axis_j(:,1);
  
  %eq. 7.207 from Computational dynamics by AA Shabana
  C=[v_i_1'*v_j;
     v_i_2'*v_j;
     v_i_1'*r_ij;
     v_i_2'*r_ij];
     
  H_i_1=[zeros(3,3) skew_sym(v_i_1)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_2=[zeros(3,3) skew_sym(v_i_2)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_j=[zeros(3,3) skew_sym(v_j)'*G_j]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_p=[eye(3) skew_sym(A_i*pt_i)'*G_i]; %eq. 7.205 from Computational dynamics by AA Shabana
  H_j_p=[eye(3) skew_sym(A_j*pt_j)'*G_j]; %eq. 7.205 from Computational dynamics by AA Shabana
  
  %eq. 7.209 from Computational dynamics by AA Shabana
  Cq=[v_j'*H_i_1, v_i_1'*H_j;...
      v_j'*H_i_2, v_i_2'*H_j;...
      r_ij'*H_i_1+v_i_1'*H_i_p, -v_i_1'*H_j_p;...
      r_ij'*H_i_2+v_i_2'*H_i_p, -v_i_2'*H_j_p];
  
  ret_joint_data.C=C;
  ret_joint_data.Cq=Cq;
  ret_joint_data.nc=4; % no. of constraints
  
  global Nc;
  global C_sys;
  global Cq_sys;
  global Cq_sys_dyn;
  global joint_index;
  global C_sys_pointer;
  global Qd_sys;
  
  joint_index=joint_index+1;
  C_sys(C_sys_pointer:C_sys_pointer+3,1)=C;
  Cq_sys(C_sys_pointer:C_sys_pointer+3,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys(C_sys_pointer:C_sys_pointer+3,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+3,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+3,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  
  q_i_dot=link_i.vel;
  q_j_dot=link_j.vel;
  q_i=link_i.cm;
  q_j=link_j.cm;

  Cqqdot_q(1,:)=[v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_1_bar,q_i,q_i_dot)+...
                (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_j*q_j_dot)'...
                (dvi_dqi(A_j,G_j,v_j_bar)'*H_i_1*q_i_dot)'+...
                 v_i_1'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot)];
         
  Cqqdot_q(2,:)=[v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_2_bar,q_i,q_i_dot)+...
                (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_j*q_j_dot)'...
                (dvi_dqi(A_j,G_j,v_j_bar)'*H_i_2*q_i_dot)'+...
                v_i_2'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot)];
         
  Cqqdot_q(3,:)=[(drij_dqi(A_i,G_i,pt_i)'*H_i_1*q_i_dot)'+...
                 r_ij'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_1_bar,q_i,q_i_dot)+...
                 (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_i_p*q_i_dot)'+...
                 v_i_1'*dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot)-...
                 (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_j_p*q_j_dot)'...
                 (drij_dqj(A_j,G_j,pt_j)'*H_i_1*q_i_dot)'-...
                 v_i_1'*dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)];
  
  Cqqdot_q(4,:)=[(drij_dqi(A_i,G_i,pt_i)'*H_i_2*q_i_dot)'+...
                  r_ij'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_2_bar,q_i,q_i_dot)+...
                 (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_i_p*q_i_dot)'+...
                  v_i_2'*dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot)-...
                 (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_j_p*q_j_dot)'...
                 (drij_dqj(A_j,G_j,pt_j)'*H_i_2*q_i_dot)'-...
                  v_i_2'*dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)];
  
  Qd_new=-Cqqdot_q*[q_i_dot;q_j_dot];
  Qd_sys(C_sys_pointer:C_sys_pointer+3,1) = Qd_new;
  
C_sys_pointer=C_sys_pointer+4;
  
  
  %Nc=nc;
endfunction

function ret_gnd_constraint = ground_constraint(ground) 
  %function to add ground constraints to the simulation
  ret_gnd_constraint.C=ground.cm;
  Cq=[eye(6) zeros(6,6)];
  ret_gnd_constraint.Cq=Cq;
  ret_gnd_constraint.nc=6; % no. of constraints 
  
  global Nc;
  global C_sys;
  global Cq_sys;
  global joint_index;
  global C_sys_pointer;
  global Qd_sys;
  global Cq_sys_dyn;
  joint_index=joint_index+1;
  C_sys(C_sys_pointer:C_sys_pointer+5,1)=ret_gnd_constraint.C;
  Cq_sys(C_sys_pointer:C_sys_pointer+5,ground.body_index*6-5:ground.body_index*6)=Cq(:,1:6);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+5,ground.body_index*6-5:ground.body_index*6)=Cq(:,1:6);
  
  Qd_sys(C_sys_pointer:C_sys_pointer+5,1) = zeros(6,1);
  C_sys_pointer=C_sys_pointer+6;
  Nc=ret_gnd_constraint.nc;

endfunction

function ret_joint_data = revolute_joint(link_i, pt_i,joint_axis_i, link_j, pt_j,joint_axis_j )
  %joint_axis_i is in BCS of body i. joint_axis_i(:,1), joint_axis_i(:,2) and joint_axis_i(:,3) form an orthogonal triad 
  %joint_axis_j is in BCS of body j
  
  %eq. 7.31 from Computational dynamics by AA Shabana
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  A_j=Ai_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  %eq. 7.72 from Computational dynamics by AA Shabana
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_j=Gi_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  
  G_i_bar=A_i'*G_i;
  G_j_bar=A_j'*G_j;
  
  r_ij=link_i.cm(1:3)+A_i*pt_i - (link_j.cm(1:3)+A_j*pt_j);
  v_i_1= A_i*joint_axis_i(:,2);
  v_i_2= A_i*joint_axis_i(:,3);
  v_j=A_j*joint_axis_j(:,1);
  v_i_1_bar= joint_axis_i(:,2);
  v_i_2_bar= joint_axis_i(:,3);
  v_j_bar=joint_axis_j(:,1);

  %eq. 7.215 from Computational dynamics by AA Shabana
  C=[r_ij;
     v_i_1'*v_j;
     v_i_2'*v_j];
     
  H_i_1=[zeros(3,3) skew_sym(v_i_1)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_2=[zeros(3,3) skew_sym(v_i_2)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_j=[zeros(3,3) skew_sym(v_j)'*G_j]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_p=[eye(3) skew_sym(A_i*pt_i)'*G_i]; %eq. 7.205 from Computational dynamics by AA Shabana
  H_j_p=[eye(3) skew_sym(A_j*pt_j)'*G_j]; %eq. 7.205 from Computational dynamics by AA Shabana
  
  %eq. 7.216 from Computational dynamics by AA Shabana
  Cq=[H_i_p, -H_j_p;...
      v_j'*H_i_1, v_i_1'*H_j;...
      v_j'*H_i_2, v_i_2'*H_j];
 
      
  ret_joint_data.C=C;
  ret_joint_data.Cq=Cq;
  ret_joint_data.nc=5; % no. of constraints
  
  global Nc;
  global C_sys;
  global Cq_sys;
  global joint_index;
  global C_sys_pointer;
  global Cq_sys_dyn;
  global Qd_sys;
  
  joint_index=joint_index+1;
  C_sys(C_sys_pointer:C_sys_pointer+4,1)=C;
  Cq_sys(C_sys_pointer:C_sys_pointer+4,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys(C_sys_pointer:C_sys_pointer+4,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+4,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+4,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  q_i_dot=link_i.vel;
  q_j_dot=link_j.vel;
  q_i=link_i.cm;
  q_j=link_j.cm;
 
  Cqqdot_q(1:3,:)=[dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot) -dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)];
  Cqqdot_q(4,1:6)=v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_1_bar,q_i,q_i_dot)+(dvi_dqi(A_i,G_i,v_i_1_bar)'*H_j*q_j_dot)';
  Cqqdot_q(4,7:12)=(dvi_dqi(A_j,G_j,v_j_bar)'*H_i_1*q_i_dot)'+v_i_1_bar'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot);
  Cqqdot_q(5,1:6)=v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_2_bar,q_i,q_i_dot)+(dvi_dqi(A_i,G_i,v_i_2_bar)'*H_j*q_j_dot)';
  Cqqdot_q(5,7:12)=(dvi_dqi(A_j,G_j,v_j_bar)'*H_i_2*q_i_dot)'+v_i_2_bar'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot);
  
  Qd=-Cqqdot_q*[q_i_dot;q_j_dot];
  Qd_sys(C_sys_pointer:C_sys_pointer+4,1) = Qd;
  C_sys_pointer=C_sys_pointer+5;
  
  %Nc=nc;
endfunction

function ret_joint_data = prismatic_joint(link_i, pt_i,joint_axis_i, link_j, pt_j,joint_axis_j )
  %joint_axis_i is in BCS of body i. joint_axis_i(:,1), joint_axis_i(:,2) and joint_axis_i(:,3) form an orthogonal triad 
  %joint_axis_j is in BCS of body j. joint_axis_j(:,2) and joint_axis_i(:,2) are perpendicular to each other
  %joint_axis_j(:,2) and joint_axis_i(:,2) are always perpendicular to each other restricting relative rotation
  
  %eq. 7.31 from Computational dynamics by AA Shabana
  A_i=Ai_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  A_j=Ai_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  %eq. 7.72 from Computational dynamics by AA Shabana
  G_i=Gi_mat(link_i.cm(4),link_i.cm(5),link_i.cm(6));
  G_j=Gi_mat(link_j.cm(4),link_j.cm(5),link_j.cm(6));
  G_i_bar=A_i'*G_i;
  G_j_bar=A_j'*G_j;
  
  r_ij=link_i.cm(1:3)+A_i*pt_i - (link_j.cm(1:3)+A_j*pt_j);
  v_i_1= A_i*joint_axis_i(:,2);
  v_i_2= A_i*joint_axis_i(:,3);
  v_j=A_j*joint_axis_j(:,1);
  h_i=A_i*joint_axis_i(:,2);
  h_j=A_j*joint_axis_j(:,3);
  
  v_i_1_bar= joint_axis_i(:,2);
  v_i_2_bar= joint_axis_i(:,3);
  v_j_bar=joint_axis_j(:,1);
  h_i_bar=joint_axis_i(:,2);
  h_j_bar=joint_axis_j(:,3);
  
  %eq. 7.207 from Computational dynamics by AA Shabana
  C=[v_i_1'*v_j;
     v_i_2'*v_j;
     v_i_1'*r_ij;
     v_i_2'*r_ij;
     h_i'*h_j];
     
  H_i_1=[zeros(3,3) skew_sym(v_i_1)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_2=[zeros(3,3) skew_sym(v_i_2)'*G_i]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_j=[zeros(3,3) skew_sym(v_j)'*G_j]; %eq. 7.210 from Computational dynamics by AA Shabana
  H_i_p=[eye(3) skew_sym(A_i*pt_i)'*G_i]; %eq. 7.205 from Computational dynamics by AA Shabana
  H_j_p=[eye(3) skew_sym(A_j*pt_j)'*G_j]; %eq. 7.205 from Computational dynamics by AA Shabana
  H_i_h=[zeros(3,3) skew_sym(h_i)'*G_i]; %eq. 7.220 from Computational dynamics by AA Shabana
  H_j_h=[zeros(3,3) skew_sym(h_j)'*G_j]; %eq. 7.220 from Computational dynamics by AA Shabana
  
  
  %eq. 7.209 from Computational dynamics by AA Shabana
  Cq=[v_j'*H_i_1, v_i_1'*H_j;...
      v_j'*H_i_2, v_i_2'*H_j;...
      r_ij'*H_i_1+v_i_1'*H_i_p, -v_i_1'*H_j_p;...
      r_ij'*H_i_2+v_i_2'*H_i_p, -v_i_2'*H_j_p;...
      h_j'*H_i_h, h_i'*H_j_h];
  
  ret_joint_data.C=C;
  ret_joint_data.Cq=Cq;
  ret_joint_data.nc=5; % no. of constraints
  
  global Nc;
  global C_sys;
  global Cq_sys;
  global Cq_sys_dyn;
  global Qd_sys;
  global joint_index;
  global C_sys_pointer;
  joint_index=joint_index+1;
  C_sys(C_sys_pointer:C_sys_pointer+4,1)=C;
  Cq_sys(C_sys_pointer:C_sys_pointer+4,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys(C_sys_pointer:C_sys_pointer+4,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+4,link_i.body_index*6-5:link_i.body_index*6)=Cq(:,1:6);
  Cq_sys_dyn(C_sys_pointer:C_sys_pointer+4,link_j.body_index*6-5:link_j.body_index*6)=Cq(:,7:12);
  
  q_i_dot=link_i.vel;
  q_j_dot=link_j.vel;
  q_i=link_i.cm;
  q_j=link_j.cm;
  Cqqdot_q(1,:)=[v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_1_bar,q_i,q_i_dot)+...
                (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_j*q_j_dot)'...
                (dvi_dqi(A_j,G_j,v_j_bar)'*H_i_1*q_i_dot)'+...
                 v_i_1'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot)];
         
  Cqqdot_q(2,:)=[v_j'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_2_bar,q_i,q_i_dot)+...
                (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_j*q_j_dot)'...
                (dvi_dqi(A_j,G_j,v_j_bar)'*H_i_2*q_i_dot)'+...
                v_i_2'*dHqidot_dqi(A_j,G_j,G_j_bar,v_j_bar,q_j,q_j_dot)];
         
  Cqqdot_q(3,:)=[(drij_dqi(A_i,G_i,pt_i)'*H_i_1*q_i_dot)'+...
                 r_ij'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_1_bar,q_i,q_i_dot)+...
                 (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_i_p*q_i_dot)'+...
                 v_i_1'*dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot)-...
                 (dvi_dqi(A_i,G_i,v_i_1_bar)'*H_j_p*q_j_dot)'...
                 (drij_dqj(A_j,G_j,pt_j)'*H_i_1*q_i_dot)'-...
                 v_i_1'*dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)];
  
  Cqqdot_q(4,:)=[(drij_dqi(A_i,G_i,pt_i)'*H_i_2*q_i_dot)'+...
                  r_ij'*dHqidot_dqi(A_i,G_i,G_i_bar,v_i_2_bar,q_i,q_i_dot)+...
                 (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_i_p*q_i_dot)'+...
                  v_i_2'*dHqidot_dqi(A_i,G_i,G_i_bar,pt_i,q_i,q_i_dot)-...
                 (dvi_dqi(A_i,G_i,v_i_2_bar)'*H_j_p*q_j_dot)'...
                 (drij_dqj(A_j,G_j,pt_j)'*H_i_2*q_i_dot)'-...
                  v_i_2'*dHqidot_dqi(A_j,G_j,G_j_bar,pt_j,q_j,q_j_dot)];
  
  Cqqdot_q(5,1:6)=h_j'*dHqidot_dqi(A_i,G_i,G_i_bar,h_i_bar,q_i,q_i_dot)+...
                 (drij_dqi(A_i,G_i,h_i_bar)'*H_j_h*q_j_dot)';
  Cqqdot_q(5,7:12)=(drij_dqi(A_j,G_j,h_j_bar)'*H_i_h*q_i_dot)'+...
                   h_i'*dHqidot_dqi(A_j,G_j,G_j_bar,h_j_bar,q_j,q_j_dot);
  
  Qd_new=-Cqqdot_q*[q_i_dot;q_j_dot];
  Qd_sys(C_sys_pointer:C_sys_pointer+4,1)=Qd_new;
  C_sys_pointer=C_sys_pointer+5;
  %Nc=nc;
endfunction

function mat_tr=row_transpose(mat,tr_info)
 mat_tr=mat;
 for i=1:size(tr_info)(2)
   if tr_info(i)==i
     continue;
   endif
   temp=mat(tr_info(i),:); % store the ith independent coordinate in temp
   j=tr_info(i);
   while j>=i+1
     mat_tr(j,:)=mat_tr(j-1,:);
     j=j-1;
   endwhile
   mat_tr(i,:)=temp;
 endfor
endfunction

function mat_tr=arow_transpose(mat,tr_info)
  mat_tr=mat;
  mat_size=size(mat)(1);
  ctr=1;
  %below looks are to determine the posistion of the dependent coordinates in the untransformed matrix
  %---------------------------
  for i=1:size(tr_info)(2)
    if i==tr_info(i)
      continue,
    endif
    for j=tr_info(i-1)+1:tr_info(i)-1
      a_tr_info(ctr)=j;
      ctr=ctr+1;
    endfor
  endfor 
  for j=tr_info(size(tr_info)(2))+1:mat_size 
    a_tr_info(ctr)=j;
    ctr=ctr+1;
  endfor
  %----------------------------
 ctr=1;
 
 for i=mat_size-size(a_tr_info)(2)+1:mat_size
   temp=mat(i,:);
   for j=i-1:-1:a_tr_info(ctr)
     mat_tr(j+1,:)=mat_tr(j,:);
   endfor
   mat_tr(a_tr_info(ctr))=temp;
   ctr=ctr+1;
 endfor
  
endfunction

function mat_tr=column_transpose(mat,tr_info)
  mat_tr=mat;
  for i=1:size(tr_info)(2)
   if tr_info(i)==i
     continue;
   endif
   temp=mat(:,tr_info(i)); % store the ith independent coordinate in temp
   j=tr_info(i);
   while j>=i+1
     mat_tr(:,j)=mat_tr(:,j-1);
     j=j-1;
   endwhile
   mat_tr(:,i)=temp;
 endfor
endfunction

function P = drem(Wmat,Vi,A_eq,mu,erest,alpha)
  %Vi = initial velocity of approach (should be negative)
  %P = impulse (to be calculated)
  %mu= coeff. of friction (1 d matrix of nc length)
  %alpha = collision parameter
  % erest = coeff. of restitution
  
  nc=length(mu); %length is nc
 % Aineq=zeros(6*nc+1+35,3*nc+35);
 % bineq=zeros(6*nc+1+35,1);
  % A*P>=0 - condition for non-tensile impulses
  A = zeros(nc,3*nc+35);
  b=zeros(nc,1);
  for i=1:nc
    A(i,3*i-2)=1;
  end
  disp("Size of A and b");
  size(A);
  size(b);
  A;
  b;
  % A*Wmat*P >= -A*V_i % Bodies should be seperating after collision
  Anp=-A*Wmat;
  bnp=zeros(nc,1);
  
  for i=1:nc
    bnp(i,1)=Vi(3*i-2);
  end
  disp("Size of Anp and bnp");
  size(Anp);
  size(bnp);
  Anp;
  bnp;
  
  % friction B*P<=0
  B=zeros(4*nc,3*nc+35);
  bf=zeros(4*nc,1);
 
  for i=1:nc
    B(i,3*i-1)=1;
    B(i,3*i-2)=-mu(i);
    B(nc+i,3*i-1)=-1;
    
    B(nc+i,3*i-2)=-mu(i);
    B(2*nc+i,3*i)=1;
    B(2*nc+i,3*i-2)=-mu(i);
    
    B(3*nc+i,3*i)=-1;
    B(3*nc+i,3*i-2)=-mu(i);
  end
  disp("Size of B and bf");
  size(B);
  size(bf);
 
  % Ad*P <= bd
  d1=zeros(1,nc);
  d2=zeros(1,3*nc+35);
  for i=1:nc
    d1(i)=1/(Wmat(3*i-2,3*i-2)^(alpha(i)));
    d2(3*i-2)=(1+erest(i))*d1(i);
  end
  Ad=-d1*A*Wmat;
  bd=d2*Vi;
  disp("Size of Ad and bd");
  size(Ad);
  size(bd);
  

  Aineq=[-A;Anp;B;Ad];
  bineq=[b;bnp;bf;bd];
  disp("Size of Aineq and bineq");
  size(Aineq);
  size(bineq);
  
  %b_eq=zeros(3*nc+35,1);
  b_eq=zeros(35,1);
  disp("Size of Aeq and beq");
  size(A_eq)
  size(b_eq)
  
  [P, obj, info, lambda] = qp (10000*rand(1,3*nc+35)',Wmat,Vi,A_eq,b_eq,[],[],[],Aineq,bineq,optimset('MaxIter',1000));
  info
endfunction

function [w,z,retcode] = LCPSolve(M,q,pivtol,maxits)
if nargin<3, pivtol = 1e-8; maxits = 1e4; end;
if nargin<4, maxits = 1e3; end;
n = length(q);
if size(M)~=[n n]; error('Matrices are not compatible'); end;
rayTerm = false;
loopcount = 0;
if min(q)>=0 % If all elements are positive a trivial solution exists
    %  As w - Mz = q, if q >= 0 then w = q a, and z = 0
    w = q;
    z = zeros(size(q));
else 
    dimen = size(M,1); % Number of rows
    % Create initial tableau
    tableau = [eye(dimen), -M, -ones(dimen, 1), q];
    % Let artificial variable enter the basis
    basis = 1:dimen; % A set of row indices in the tableau
    [~,locat] = min(tableau(:,end)); % Row of minimum element in column 
                                     % 2*dimen+1 (last of tableau)
    basis(locat) = 2*dimen+1; % Replace that index with the column
    cand = locat + dimen;
    pivot = tableau(locat,:)/tableau(locat,2*dimen+1);
    % From each column subtract the column 2*dimen+1, multiplied by pivot
    tableau = tableau - tableau(:,2*dimen+1)*pivot; 
    tableau(locat,:) = pivot; % set all elements of row locat to pivot
    % Perform complementary pivoting
    while max(basis) == 2*dimen+1 && loopcount < maxits
        loopcount = loopcount + 1;
        eMs = tableau(:,cand); % This is used to check convergence (pivtol)
        missmask = eMs <= 0;  % Check if elements of eMs are less than zero
        quots = tableau(:,2*dimen+2)./eMs;
        quots(missmask) = Inf;
        [~,locat] = min(quots);
        % Check if at least one element is not missing
        if  sum(missmask)~=dimen && abs(eMs(locat)) > pivtol 
            % Reduce tableau
            pivot = tableau(locat,:)/tableau(locat,cand);
            tableau = tableau - tableau(:,cand)*pivot;
            tableau(locat,:) = pivot;
            oldVar = basis(locat);
            % New variable enters the basis
            basis(locat) = cand;
            % Select next candidate for entering the basis
            if oldVar > dimen
                cand = oldVar - dimen;
            else
                cand = oldVar + dimen;
            end
        else
            rayTerm = true; % Problem was solved
            break % Break out of the while loop
        end
    end
    % Return the solution to the LCP
    vars = zeros(2*dimen+1,1);
    vars(basis) = tableau(:,2*dimen+2).';
    w = vars(1:dimen,1);
    z = vars(dimen+1:2*dimen,1);
end
if rayTerm
    retcode = [2, loopcount];  % Ray termination
else
    retcode = [1, loopcount];  % Success
end
endfunction

function reduced_angle = ret_angle(angle)
  for i=1:size(angle)(1)
    if angle(i) <=0
      reduced_angle(i,1)=abs(angle(i))/(2*pi());
      reduced_angle(i,1)=reduced_angle(i,1)-floor(reduced_angle(i,1));
      reduced_angle(i,1)=-reduced_angle(i,1)*2*pi();
    else
      reduced_angle(i,1)=abs(angle(i,1))/(2*pi());
      reduced_angle(i,1)=reduced_angle(i,1)-floor(reduced_angle(i,1));
      reduced_angle(i,1)=reduced_angle(i,1)*2*pi();
  endif
endfor
endfunction

global Nb=12; %no. of bodies
%global Nc=6*Nb-13; % Frame, boom ,rocker, stick, tophammer, piston
global Np=24; % no. of collision points
global C_sys=zeros(6*Nb,1);
global Cq_sys_dyn; % Jacobian matrix (non square) to be used for dynamic calculation.
global Cq_sys=zeros(6*Nb,6*Nb);
global joint_index=0;
global C_sys_pointer=1;
global Qd_sys=zeros(1,1);
global Qe_sys=zeros(6*Nb,1);
global Qv_sys=zeros(6*Nb,1);


init_rot=deg2rad(30);
%define ground body

gnd_T1.mass =1E8; % mass of link
gnd_T1.I_theta_bar = [1E8 0 0; 0 1E8 0; 0 0 1E8]; 
% Coordinates of points in BCS of ground.
gnd_T1.T1=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T1.T2=Rot_y(init_rot)*[-0.5;0;1.3];
gnd_T1.T3=Rot_y(init_rot)*[-0.5;0;-1.3];
gnd_T1.T4=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T1.cm=[0;0;0;0;0;0]; %initial position in GCS
gnd_T1.vel=zeros(6,1); %initial vel in GCS
init_cond.gnd_T1=gnd_T1.cm;
gnd_T1.body_index=1;

gnd_T2.mass =1E8; % mass of link
gnd_T2.I_theta_bar = [1E8 0 0; 0 1E8 0; 0 0 1E8]; 
% Coordinates of points in BCS of ground.
gnd_T2.T1=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T2.T2=Rot_y(init_rot)*[-0.5;0;1.3];
gnd_T2.T3=Rot_y(init_rot)*[-0.5;0;-1.3];
gnd_T2.T4=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T2.cm=[0;0;0;0;0;0]; %initial position in GCS
gnd_T2.vel=zeros(6,1); %initial vel in GCS
init_cond.gnd_T2=gnd_T2.cm;
gnd_T2.body_index=2;

gnd_T3.mass =1E8; % mass of link
gnd_T3.I_theta_bar = [1E8 0 0; 0 1E8 0; 0 0 1E8]; 
% Coordinates of points in BCS of ground.
gnd_T3.T1=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T3.T2=Rot_y(init_rot)*[-0.5;0;1.3];
gnd_T3.T3=Rot_y(init_rot)*[-0.5;0;-1.3];
gnd_T3.T4=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T3.cm=[0;0;0;0;0;0]; %initial position in GCS
gnd_T3.vel=zeros(6,1); %initial vel in GCS
init_cond.gnd_T3=gnd_T3.cm;
gnd_T3.body_index=3;

gnd_T4.mass =1E8; % mass of link
gnd_T4.I_theta_bar = [1E8 0 0; 0 1E8 0; 0 0 1E8]; 
% Coordinates of points in BCS of ground.
gnd_T4.T1=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T4.T2=Rot_y(init_rot)*[-0.5;0;1.3];
gnd_T4.T3=Rot_y(init_rot)*[-0.5;0;-1.3];
gnd_T4.T4=Rot_y(init_rot)*[2.7;0;1.3];
gnd_T4.cm=[0;0;0;0;0;0]; %initial position in GCS
gnd_T4.vel=zeros(6,1); %initial vel in GCS
init_cond.gnd_T4=gnd_T4.cm;
gnd_T4.body_index=4;



% define  frame
frame.mass =15000; % mass of link

frame.I_theta_bar=[1.50E+04	5.25E+03	1.99E+02;...
                   5.25E+03	3.10E+04	1.27E+02;...
                   1.99E+02	1.27E+02	2.73E+04]; 
frame.T1=[2.7;-1.6;1.3];
frame.T2=[-0.5;-1.6;1.3];
frame.T3=[-0.5;-1.6;-1.3];
frame.T4=[2.7;-1.6;-1.3];
frame.pt1=[1.35;0.08;0.55]; %[1.35;-0.55;0.65];
frame.pt2=[1.35;-1.5;0.55];
frame.cm_BCS=[0;0;0]; % cm position in bcs

frame.cm=[-2.833404
          -2.495431
          -1.342040
          -0.075566
          -0.433024
          -0.352025]; %initial position in GCS

frame.vel=zeros(6,1); %initial vel in GCS
frame.M_mat=calc_M_i(frame);
init_cond.frame=frame.cm;
frame.body_index=5;

% define  boom
boom.mass =1000; % mass of link

boom.I_theta_bar=[1.65E+01	4.14E-01	0.00E+00;...
                   4.14E-01	1.07E+03	0.00E+00;...
                   0.00E+00	0.00E+00	1.07E+03];                   
boom.pt1=[-1.5;0;0];
boom.pt3=[0.15;-0.32;0];
boom.pt4=[0.3;0.32;0];
boom.pt5=[1.5;0;0];
boom.cm_BCS=[0;0;0];
boom.cm=[4.6827
   4.7503
   1.3568
   2.1508
   5.4459
   6.4889]+0*rand(6,1); %initial position in GCS

boom.vel=zeros(6,1); %initial vel in GCS
boom.body_index=6;
init_cond.boom=boom.cm;

% define rocker
rocker.mass=521;
rocker.I_theta_bar=[3.99E+01	-3.89E+00	-5.01E-01;...
                   -3.89E+00	5.02E+01	2.27E+00;...
                   -5.01E-01	2.27E+00	6.18E+01];

rocker.pt5=[-0.2;0;0];
rocker.pt6=[0.8/3;0.4;0];
rocker.pt7=[0;-0.2;0.1];
rocker.pt8=[0.25;-0.15;0.3];
rocker.cm_BCS=[0;0;0];
rocker.cm=ones(6,1)+randn(6,1); %initial position in GCS      
       
rocker.cm=[-0.33820
           -0.40392
            1.77061
            1.57955
            1.97426
            2.14515]+0*rand(6,1);
            
rocker.vel=zeros(6,1); %initial vel in GCS
rocker.body_index=7;
init_cond.rocker=rocker.cm;

% define stick
stick.mass=900;

stick.I_theta_bar=[2.70E+03	-2.05E+00	0.00E+00;...
                  -2.05E+00	3.13E+01	5.17E-01;...
                   0.00E+00	5.17E-01	2.69E+03];

stick.pt8=[0;0;-0.4];
stick.pt9=[0;-9.2/2;0.1];
stick.pt10=[0.25;2.95;0];
stick.cm_BCS=[0;0;0];
stick.cm=ones(6,1)+randn(6,1); %initial position in GCS

stick.cm=[2.90569
   1.99819
   2.32726
   0.21119
   0.86047
   0.24365];
       
stick.vel=zeros(6,1); %initial vel in GCS
stick.body_index=8;
init_cond.stick=stick.cm;

% define top hammer
top_hammer.mass=150;

top_hammer.I_theta_bar=[1.44E+01	-4.14E-03	0.00E+00;...
                       -4.14E-03	5.58E+00	0.00E+00;...
                        0.00E+00	0.00E+00	1.66E+01];


top_hammer.pt10=[-0.2;0.5;0];
top_hammer.pt11=[-0.07;0.2;0];
top_hammer.pt12=[-0.07;-0.181;0];
top_hammer.cp_TH_drill_rod_imp=[-0.07;-0.100;0]; % Bit in impact position
top_hammer.cp_TH_drill_rod_ext=[-0.07;-0.156;0]; % Bit in extended position
top_hammer.cm_BCS=[0;0;0];
top_hammer.cm=ones(6,1)+7*randn(6,1); %initial position in GCS
top_hammer.cm=[4.499 %5.15583
               3.30867
              -0.949714 %0.80054
               -5.80056
               1.02204
               0];%-3.12499];
     
top_hammer.vel=zeros(6,1); %initial vel in GCS
top_hammer.body_index=9;
init_cond.top_hammer=top_hammer.cm;

% define hammer piston
piston.mass=19.2;

piston.I_theta_bar=[2.61E-01	0.00E+00	0.00E+00;...
                    0.00E+00	2.06E-02	0.00E+00;...
                    0.00E+00	0.00E+00	2.61E-01];

piston.pt11=[0;0;0];
piston.cm_BCS=[0;0;0];
piston.cp_piston_tool=[0;-0.2;0]; %contact point between piston and tool
piston.cm=ones(6,1)+50*randn(6,1); %initial position in GCS
piston.cm=[11.1141
        9.3567
       19.5794
       19.5782
       68.1960
      -23.8704]+6*rand(6,1);

piston.cm(1:6)=[12.050
                12.270
                22.794
                22.355
                70.126
               -21.640];   
               
piston.cm
piston.vel=zeros(6,1); %initial vel in GCS
piston.body_index=10;
init_cond.piston=piston.cm;

% define rotary_gear
rotary_gear.mass=8.7;

rotary_gear.I_theta_bar=[2.46E-02	0.00E+00	0.00E+00;...
                         0.00E+00	4.56E-02	0.00E+00;...
                         0.00E+00	0.00E+00	2.46E-02];

rotary_gear.pt12=[0;0;0]; %joint with hammer
rotary_gear.pt13=[0;1;0]; %joint with drill rod
rotary_gear.cm_BCS=[0;0;0];
rotary_gear.cm=ones(6,1)+randn(6,1); %initial position in GCS
rotary_gear.cm=-[-1.21885
       2.08060
       0.74704
       0.74301
       0.30571
      -0.15642];    
rotary_gear.cm=[4.43074
                5.12836
                -0.92838
                -4.88561
                 0.53015
                -1.42097];   
rotary_gear.vel=zeros(6,1); %initial vel in GCS
rotary_gear.body_index=11;
init_cond.rotary_gear=rotary_gear.cm;

% define drill_rod_bit
drill_rod_bit.mass=101.5;
%{
drill_rod_bit.I_theta_bar=[200,0,0;0,500,0;0,0,100];
%}

drill_rod_bit.I_theta_bar=[2.79E+02	2.25E-06	0.00E+00;...
                           2.25E-06	6.99E-02	2.72E-05;...
                           0.00E+00	2.72E-05	2.79E+02];


drill_rod_bit.pt13=[0;5*9.2/12;0];
drill_rod_bit.pt14=[0;-5*9.2/12;0];
drill_rod_bit.cp_tool_piston=[0;3.0;0];
drill_rod_bit.cp_tool_top_hammer(:,1)=[0;2.950;0];
drill_rod_bit.cp_tool_top_hammer(:,2)=[0;2.920;0];

drill_rod_bit.cp_gauge_button(:,1)=[0.0481436;-1.99750;-0.0277957];
drill_rod_bit.cp_gauge_button(:,2)=[0.0;-1.99750;-0.0555915];
drill_rod_bit.cp_gauge_button(:,3)=[-0.0481436;-1.99750;-0.0277957];
drill_rod_bit.cp_gauge_button(:,4)=[-0.0481436;-1.99750;0.0277957];
drill_rod_bit.cp_gauge_button(:,5)=[0.0;-1.99750;0.0555915];
drill_rod_bit.cp_gauge_button(:,6)=[-0.0481436;-1.99750;0.0277957];
drill_rod_bit.cp_intermediate_button(:,1)=[0.0210479;-2.00359;-0.0364561];
drill_rod_bit.cp_intermediate_button(:,2)=[-0.0210479;-2.00359;-0.0364561];
drill_rod_bit.cp_intermediate_button(:,3)=[-0.0210479;-2.00359;0.0364561];
drill_rod_bit.cp_intermediate_button(:,4)=[0.0210479;-2.00359;0.0364561];
drill_rod_bit.cp_face_button(:,1)=[0.00792999;-2.006;-0.00792999];
drill_rod_bit.cp_face_button(:,2)=[0;-2.006;-0.028];
drill_rod_bit.cp_face_button(:,3)=[-0.0075;-2.006;0.01299904];
drill_rod_bit.cp_face_button(:,4)=[0;-2.006;0.028];
drill_rod_bit.cm_BCS=[0;0;0];
drill_rod_bit.cm=9*ones(6,1)+300*randn(6,1); %initial position in GCS
drill_rod_bit.cm=-[195.898
                  -76.723
                   495.334
                   288.432
                  -257.052
                  -468.975];
%drill_rod_bit.cm=ones(6,1)+0*randn(6,1);
drill_rod_bit.vel=zeros(6,1); %initial vel in GCS
drill_rod_bit.body_index=12;
init_cond.drill_rod_bit=drill_rod_bit.cm;


axis_frame_boom=[0 0 1;... % joint axis and orthogonal triad in BCS of frame
                 0 1 0;...
                 1 0 0];
axis_boom_frame=[0 0 1;... % joint axis and orthogonal triad in BCS of boom
                 0 1 0;...
                 1 0 0];
axis_boom_rocker=[0 0 1;... % joint axis and orthogonal triad in BCS of boom
                  0 1 0;...
                  1 0 0];
axis_rocker_boom=[0 0 1;... % joint axis and orthogonal triad in BCS of rocker
                  0 1 0;...
                  1 0 0];
                  
axis_rocker_stick=[1 0 0;... % joint axis and orthogonal triad in BCS of rocker
                   0 1 0;...
                   0 0 1];
axis_stick_rocker=[1 0 0;... % joint axis and orthogonal triad in BCS of stick
                   0 1 0;...
                   0 0 1];
                   
axis_stick_hammer=[0 1 0;... % joint axis and orthogonal triad in BCS of rocker
                   1 0 0;...
                   0 0 1];
axis_hammer_stick=[0 1 0;... % joint axis and orthogonal triad in BCS of stick
                   1 0 0;...
                   0 0 1];
                   
axis_hammer_piston=[0 0 1;... % joint axis and orthogonal triad in BCS of hammer
                    1 0 0;...
                    0 1 0];
axis_piston_hammer=[0 0 1;... % joint axis and orthogonal triad in BCS of piston
                    1 0 0;...
                    0 1 0];
                  
axis_hammer_gear=[0 0 1;... % joint axis and orthogonal triad in BCS of rocker
                  1 0 0;...
                  0 1 0];
axis_gear_hammer=[0 1 0;... % joint axis and orthogonal triad in BCS of stick
                  1 0 0;...
                  0 0 1];
                      
axis_gear_drillrod=[0 0 1;... % joint axis and orthogonal triad in BCS of rocker
                    1 0 0;...
                    0 1 0];
axis_drillrod_gear=[0 0 1;... % joint axis and orthogonal triad in BCS of stick
                    1 0 0;...
                    0 1 0];
                  
                  
A_ref=Rot_y(init_rot);
theta_ref=acos(A_ref(3,3));
phi_ref=acos(-A_ref(2,3)/sin(theta_ref));
psi_ref=-acos(A_ref(3,2)/sin(theta_ref));
A_ref_2=Ai_mat(phi_ref,theta_ref,psi_ref);

A_ref_boom=Rot_y(init_rot)*Rot_z(deg2rad(45));
theta_ref_boom=acos(A_ref_boom(3,3));
phi_ref_boom=acos(-A_ref_boom(2,3)/sin(theta_ref_boom));
psi_ref_boom=-acos(A_ref_boom(3,2)/sin(theta_ref_boom));
A_ref_boom_2=Ai_mat(phi_ref_boom,theta_ref_boom,psi_ref_boom);

A_ref_stick=Rot_y(init_rot)*Rot_x(deg2rad(5));
theta_ref_stick=acos(A_ref_stick(3,3));
phi_ref_stick=acos(-A_ref_stick(2,3)/sin(theta_ref_stick));
psi_ref_stick=-acos(A_ref_stick(3,2)/sin(theta_ref_stick));
A_ref_stick_2=Ai_mat(phi_ref_stick,theta_ref_stick,psi_ref_stick);

A_ref_rotary_gear=Rot_y(init_rot)*Rot_x(deg2rad(5))*Rot_y(deg2rad(0));
theta_ref_rotary_gear=acos(A_ref_rotary_gear(3,3));
phi_ref_rotary_gear=acos(-A_ref_rotary_gear(2,3)/sin(theta_ref_rotary_gear));
psi_ref_rotary_gear=-acos(A_ref_rotary_gear(3,2)/sin(theta_ref_rotary_gear));

ctr=1;
loop_ctr=1;

for ctr=1:1000
del_q=zeros(6*Nb,1);
C_sys_pointer=1;
C_sys=zeros(6*Nb,6*Nb)*C_sys;
Cq_sys=zeros(6*Nb,6*Nb)*Cq_sys;

revolute_joint(frame,frame.pt1,axis_frame_boom,boom,boom.pt1,axis_boom_frame);
revolute_joint(boom,boom.pt5,axis_boom_rocker,rocker,rocker.pt5,axis_rocker_boom);
revolute_joint(rocker,rocker.pt8,axis_rocker_stick,stick,stick.pt8,axis_stick_rocker);
prismatic_joint(stick,stick.pt10,axis_stick_hammer,top_hammer,top_hammer.pt10,axis_hammer_stick);
prismatic_joint(top_hammer,top_hammer.pt11,axis_hammer_piston,piston,piston.pt11,axis_piston_hammer);
revolute_joint(top_hammer,top_hammer.pt12,axis_hammer_gear,rotary_gear,rotary_gear.pt12,axis_gear_hammer);
prismatic_joint(rotary_gear,rotary_gear.pt13,axis_gear_drillrod,drill_rod_bit,drill_rod_bit.pt13,axis_drillrod_gear);

Cq_sys;

% initial cond. for gnd_T4
link_i=gnd_T1;     
C_sys(C_sys_pointer:C_sys_pointer+2)=link_i.cm(1:3)-[0;0;0];
C_sys(C_sys_pointer+3)=link_i.cm(4)-0;
C_sys(C_sys_pointer+4)=link_i.cm(5)-0;
C_sys(C_sys_pointer+5)=link_i.cm(6)-0;
Cq_sys(C_sys_pointer:C_sys_pointer+5,link_i.body_index*6-5:link_i.body_index*6)=eye(6,6);
C_sys_pointer=C_sys_pointer+6;
Cq_sys;

% initial cond. for gnd_T4
link_i=gnd_T2;     
C_sys(C_sys_pointer:C_sys_pointer+2)=link_i.cm(1:3)-[0;0;0];
C_sys(C_sys_pointer+3)=link_i.cm(4)-0;
C_sys(C_sys_pointer+4)=link_i.cm(5)-0;
C_sys(C_sys_pointer+5)=link_i.cm(6)-0;
Cq_sys(C_sys_pointer:C_sys_pointer+5,link_i.body_index*6-5:link_i.body_index*6)=eye(6,6);
C_sys_pointer=C_sys_pointer+6;
Cq_sys;

% initial cond. for gnd_T4
link_i=gnd_T3;     
C_sys(C_sys_pointer:C_sys_pointer+2)=link_i.cm(1:3)-[0;0;0];
C_sys(C_sys_pointer+3)=link_i.cm(4)-0;
C_sys(C_sys_pointer+4)=link_i.cm(5)-0;
C_sys(C_sys_pointer+5)=link_i.cm(6)-0;
Cq_sys(C_sys_pointer:C_sys_pointer+5,link_i.body_index*6-5:link_i.body_index*6)=eye(6,6);
C_sys_pointer=C_sys_pointer+6;
Cq_sys;

% initial cond. for gnd_T4
link_i=gnd_T4;     
C_sys(C_sys_pointer:C_sys_pointer+2)=link_i.cm(1:3)-[0;0;0];
C_sys(C_sys_pointer+3)=link_i.cm(4)-0;
C_sys(C_sys_pointer+4)=link_i.cm(5)-0;
C_sys(C_sys_pointer+5)=link_i.cm(6)-0;
Cq_sys(C_sys_pointer:C_sys_pointer+5,link_i.body_index*6-5:link_i.body_index*6)=eye(6,6);
C_sys_pointer=C_sys_pointer+6;
Cq_sys;

% initial cond. for frame
link_i=frame;     
C_sys(C_sys_pointer:C_sys_pointer+2)=link_i.cm(1:3)-[0;1.25;0];
C_sys(C_sys_pointer+3)=link_i.cm(4)-phi_ref;
C_sys(C_sys_pointer+4)=link_i.cm(5)-theta_ref;
C_sys(C_sys_pointer+5)=link_i.cm(6)-psi_ref;
Cq_sys(C_sys_pointer:C_sys_pointer+5,link_i.body_index*6-5:link_i.body_index*6)=eye(6,6);
C_sys_pointer=C_sys_pointer+6;
Cq_sys;

% initial cond. for boom
link_i=boom;
C_sys(C_sys_pointer,1)=link_i.cm(6)-deg2rad(-45);
Cq_sys(C_sys_pointer,link_i.body_index*6)=1;

%C_sys(C_sys_pointer,1)=link_i.cm(2)-0.25;
%Cq_sys(C_sys_pointer,link_i.body_index*6-4)=1;

C_sys_pointer=C_sys_pointer+1;
Cq_sys;

% initial cond. for rocker
link_i=rocker;
C_sys(C_sys_pointer,1)=link_i.cm(6)-deg2rad(90);
Cq_sys(C_sys_pointer,link_i.body_index*6)=1;
C_sys_pointer=C_sys_pointer+1;


% initial cond. for stick
link_i=stick;
link_j=rocker;
C_sys(C_sys_pointer,1)=link_i.cm(6)-psi_ref_stick;
Cq_sys(C_sys_pointer,link_i.body_index*6)=1;
C_sys_pointer=C_sys_pointer+1;


% initial cond. for hammer
link_i=top_hammer;
C_sys(C_sys_pointer,1)=link_i.cm(2)-5.30867;
Cq_sys(C_sys_pointer,link_i.body_index*6-4)=1;
C_sys_pointer=C_sys_pointer+1;



% initial cond. for piston
link_i=piston;
C_sys(C_sys_pointer,1)=link_i.cm(2)-5.45809; %5.52771;
Cq_sys(C_sys_pointer,link_i.body_index*6-4)=1;
C_sys_pointer=C_sys_pointer+1;

link_i=rotary_gear;
C_sys(C_sys_pointer,1)=link_i.cm(4)-phi_ref_rotary_gear;
Cq_sys(C_sys_pointer,link_i.body_index*6-2)=1;

C_sys_pointer=C_sys_pointer+1;

% initial cond. for drill_rod_bit
link_i=drill_rod_bit;
C_sys(C_sys_pointer,1)=link_i.cm(2)-2.27026; %2.26927;
Cq_sys(C_sys_pointer,link_i.body_index*6-4)=1;
C_sys_pointer=C_sys_pointer+1;
C_sys;

del_q=-inv(Cq_sys)*C_sys;
C_sys_prev=C_sys;

if isnan(del_q) 
  gnd_T1.cm=init_cond.gnd_T1;
  gnd_T2.cm=init_cond.gnd_T2;
  gnd_T3.cm=init_cond.gnd_T3;
  gnd_T4.cm=init_cond.gnd_T4;
  frame.cm=init_cond.frame;
  boom.cm=init_cond.boom;
  rocker.cm=init_cond.rocker;
  stick.cm=init_cond.stick+2*randn(6,1);
  init_cond.stick=stick.cm;
  top_hammer.cm=init_cond.top_hammer;
  piston.cm=init_cond.piston;
  rotary_gear.cm=init_cond.rotary_gear;
  drill_rod_bit.cm=init_cond.drill_rod_bit;
  ctr=1;
  loop_ctr=loop_ctr+1;
  if loop_ctr>=100
    disp("Did not converge due to Nan values");
    break;
  endif
  continue;
endif
gnd_T1.cm=gnd_T1.cm+del_q(gnd_T1.body_index*6-5:gnd_T1.body_index*6);
gnd_T2.cm=gnd_T2.cm+del_q(gnd_T2.body_index*6-5:gnd_T2.body_index*6);
gnd_T3.cm=gnd_T3.cm+del_q(gnd_T3.body_index*6-5:gnd_T3.body_index*6);
gnd_T4.cm=gnd_T4.cm+del_q(gnd_T4.body_index*6-5:gnd_T4.body_index*6);
frame.cm=frame.cm+del_q(frame.body_index*6-5:frame.body_index*6);
boom.cm=boom.cm+del_q(boom.body_index*6-5:boom.body_index*6);
rocker.cm=rocker.cm+del_q(rocker.body_index*6-5:rocker.body_index*6);
stick.cm=stick.cm+del_q(stick.body_index*6-5:stick.body_index*6);
top_hammer.cm=top_hammer.cm+del_q(top_hammer.body_index*6-5:top_hammer.body_index*6);
piston.cm=piston.cm+del_q(piston.body_index*6-5:piston.body_index*6);
rotary_gear.cm=rotary_gear.cm+del_q(rotary_gear.body_index*6-5:rotary_gear.body_index*6);
drill_rod_bit.cm=drill_rod_bit.cm+del_q(drill_rod_bit.body_index*6-5:drill_rod_bit.body_index*6);

if rem(ctr,50) == 0
norm(del_q)
endif 

q(1:6,1)=gnd_T1.cm;
q(7:12,1)=gnd_T2.cm;
q(13:18,1)=gnd_T3.cm;
q(19:24,1)=gnd_T4.cm;
q(25:30,1)=frame.cm;
q(31:36,1)=boom.cm;
q(37:42,1)=rocker.cm;
q(43:48,1)=stick.cm;
q(49:54,1)=top_hammer.cm;
q(55:60,1)=piston.cm;
q(61:66,1)=rotary_gear.cm;
q(67:72,1)=drill_rod_bit.cm;
convergence_mat(ctr)=norm(del_q);
ctr=ctr+1;
if norm(del_q)<1e-5
norm(del_q)  
disp("Converged");
break;
endif
endfor
C_sys;
q;
plot(log(convergence_mat));

% Construct M matrix
M_sys(6*gnd_T1.body_index-5:6*gnd_T1.body_index-3,6*gnd_T1.body_index-5:6*gnd_T1.body_index-3)=gnd_T1.mass*eye(3,3);
M_sys(6*gnd_T1.body_index-2:6*gnd_T1.body_index,6*gnd_T1.body_index-2:6*gnd_T1.body_index)=gnd_T1.I_theta_bar; 
M_sys(6*gnd_T2.body_index-5:6*gnd_T2.body_index-3,6*gnd_T2.body_index-5:6*gnd_T2.body_index-3)=gnd_T2.mass*eye(3,3);
M_sys(6*gnd_T2.body_index-2:6*gnd_T2.body_index,6*gnd_T2.body_index-2:6*gnd_T2.body_index)=gnd_T2.I_theta_bar; 
M_sys(6*gnd_T3.body_index-5:6*gnd_T3.body_index-3,6*gnd_T3.body_index-5:6*gnd_T3.body_index-3)=gnd_T3.mass*eye(3,3);
M_sys(6*gnd_T3.body_index-2:6*gnd_T3.body_index,6*gnd_T3.body_index-2:6*gnd_T3.body_index)=gnd_T3.I_theta_bar; 
M_sys(6*gnd_T4.body_index-5:6*gnd_T4.body_index-3,6*gnd_T4.body_index-5:6*gnd_T4.body_index-3)=gnd_T4.mass*eye(3,3);
M_sys(6*gnd_T4.body_index-2:6*gnd_T4.body_index,6*gnd_T4.body_index-2:6*gnd_T4.body_index)=gnd_T4.I_theta_bar; 
M_sys(6*frame.body_index-5:6*frame.body_index-3,6*frame.body_index-5:6*frame.body_index-3)=frame.mass*eye(3,3);
M_sys(6*frame.body_index-2:6*frame.body_index,6*frame.body_index-2:6*frame.body_index)=frame.I_theta_bar; 
M_sys(6*boom.body_index-5:6*boom.body_index-3,6*boom.body_index-5:6*boom.body_index-3)=boom.mass*eye(3,3);
M_sys(6*boom.body_index-2:6*boom.body_index,6*boom.body_index-2:6*boom.body_index)=boom.I_theta_bar; 
M_sys(6*rocker.body_index-5:6*rocker.body_index-3,6*rocker.body_index-5:6*rocker.body_index-3)=rocker.mass*eye(3,3);
M_sys(6*rocker.body_index-2:6*rocker.body_index,6*rocker.body_index-2:6*rocker.body_index)=rocker.I_theta_bar;
M_sys(6*stick.body_index-5:6*stick.body_index-3,6*stick.body_index-5:6*stick.body_index-3)=stick.mass*eye(3,3);
M_sys(6*stick.body_index-2:6*stick.body_index,6*stick.body_index-2:6*stick.body_index)=stick.I_theta_bar;
M_sys(6*top_hammer.body_index-5:6*top_hammer.body_index-3,6*top_hammer.body_index-5:6*top_hammer.body_index-3)=top_hammer.mass*eye(3,3);
M_sys(6*top_hammer.body_index-2:6*top_hammer.body_index,6*top_hammer.body_index-2:6*top_hammer.body_index)=top_hammer.I_theta_bar;
M_sys(6*piston.body_index-5:6*piston.body_index-3,6*piston.body_index-5:6*piston.body_index-3)=piston.mass*eye(3,3);
M_sys(6*piston.body_index-2:6*piston.body_index,6*piston.body_index-2:6*piston.body_index)=piston.I_theta_bar;
M_sys(6*rotary_gear.body_index-5:6*rotary_gear.body_index-3,6*rotary_gear.body_index-5:6*rotary_gear.body_index-3)=rotary_gear.mass*eye(3,3);
M_sys(6*rotary_gear.body_index-2:6*rotary_gear.body_index,6*rotary_gear.body_index-2:6*rotary_gear.body_index)=rotary_gear.I_theta_bar;
M_sys(6*drill_rod_bit.body_index-5:6*drill_rod_bit.body_index-3,6*drill_rod_bit.body_index-5:6*drill_rod_bit.body_index-3)=drill_rod_bit.mass*eye(3,3);
M_sys(6*drill_rod_bit.body_index-2:6*drill_rod_bit.body_index,6*drill_rod_bit.body_index-2:6*drill_rod_bit.body_index)=drill_rod_bit.I_theta_bar;


% Construct J_mbd matrix 
body_contact_no=1; % This counts each contact for each body. So if there are two bodies contacting at 1 point this variable will be 2 in total
system_contact_ctr=1;  % contact the no. of contact points in the system

%gnd_T1
Gi_body=Gi_mat(gnd_T1.cm(4),gnd_T1.cm(5),gnd_T1.cm(6));
Ai_body=Ai_mat(gnd_T1.cm(4),gnd_T1.cm(5),gnd_T1.cm(6));
body_index=gnd_T1.body_index;
u_i_bar_body=gnd_T1.T1;
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=-eye(3,3); 
body_contact_no=body_contact_no+3;
gnd_T1_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;

%gnd_T2
Gi_body=Gi_mat(gnd_T2.cm(4),gnd_T2.cm(5),gnd_T2.cm(6));
Ai_body=Ai_mat(gnd_T2.cm(4),gnd_T2.cm(5),gnd_T2.cm(6));
body_index=gnd_T2.body_index;
u_i_bar_body=gnd_T2.T2;
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=-eye(3,3); 
body_contact_no=body_contact_no+3;
gnd_T2_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;

%gnd_T3
Gi_body=Gi_mat(gnd_T3.cm(4),gnd_T3.cm(5),gnd_T3.cm(6));
Ai_body=Ai_mat(gnd_T3.cm(4),gnd_T3.cm(5),gnd_T3.cm(6));
body_index=gnd_T3.body_index;
u_i_bar_body=gnd_T3.T3;
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=-eye(3,3); 
body_contact_no=body_contact_no+3;
gnd_T3_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;

%gnd_T4
Gi_body=Gi_mat(gnd_T4.cm(4),gnd_T4.cm(5),gnd_T4.cm(6));
Ai_body=Ai_mat(gnd_T4.cm(4),gnd_T4.cm(5),gnd_T4.cm(6));
body_index=gnd_T4.body_index;
u_i_bar_body=gnd_T4.T4;
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=-eye(3,3); 
body_contact_no=body_contact_no+3;
gnd_T4_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;
disp("gnd_T4");
disp((gnd_T4.cm(1:3)+Ai_body*u_i_bar_body)');


%frame points
frame.T1=[2.7;-1.25;1.3];%[2.7;-1.6;1.3];
frame.T2=[-0.5;-1.25;1.3];%[-0.5;-1.6;0];
frame.T3=[-0.5;-1.25;-1.3];%[2.7;-1.6;-1.3];
frame.T4=[2.7;-1.25;-1.3];
Gi_body=Gi_mat(frame.cm(4),frame.cm(5),frame.cm(6));
Ai_body=Ai_mat(frame.cm(4),frame.cm(5),frame.cm(6));

J_mbd(6*frame.body_index-5:6*frame.body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*frame.body_index-2:6*frame.body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*frame.T1)';
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-5:6*frame.body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-2:6*frame.body_index)=-skew_sym(Ai_body*frame.T1)*Gi_body;
R_mbd(gnd_T1_contact_no:gnd_T1_contact_no+2,gnd_T1_contact_no:gnd_T1_contact_no+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,gnd_T1_contact_no:gnd_T1_contact_no+2)=eye(3,3);
body_contact_no=body_contact_no+3;
system_contact_ctr=system_contact_ctr;
disp("Frame T1: ");
disp((frame.cm(1:3)+Ai_body*frame.T1)');


J_mbd(6*frame.body_index-5:6*frame.body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*frame.body_index-2:6*frame.body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*frame.T2)';
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-5:6*frame.body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-2:6*frame.body_index)=-skew_sym(Ai_body*frame.T2)*Gi_body;
R_mbd(gnd_T2_contact_no:gnd_T2_contact_no+2,gnd_T2_contact_no:gnd_T2_contact_no+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,gnd_T2_contact_no:gnd_T2_contact_no+2)=eye(3,3);
body_contact_no=body_contact_no+3;
system_contact_ctr=system_contact_ctr;
disp("Frame T2: ");
disp((frame.cm(1:3)+Ai_body*frame.T2)');

J_mbd(6*frame.body_index-5:6*frame.body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*frame.body_index-2:6*frame.body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*frame.T3)';
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-5:6*frame.body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-2:6*frame.body_index)=-skew_sym(Ai_body*frame.T3)*Gi_body;
R_mbd(gnd_T3_contact_no:gnd_T3_contact_no+2,gnd_T3_contact_no:gnd_T3_contact_no+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,gnd_T3_contact_no:gnd_T3_contact_no+2)=eye(3,3);
body_contact_no=body_contact_no+3;
system_contact_ctr=system_contact_ctr;
disp("Frame T3: ");
disp((frame.cm(1:3)+Ai_body*frame.T3)');


J_mbd(6*frame.body_index-5:6*frame.body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*frame.body_index-2:6*frame.body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*frame.T4)';
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-5:6*frame.body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*frame.body_index-2:6*frame.body_index)=-skew_sym(Ai_body*frame.T4)*Gi_body;
R_mbd(gnd_T4_contact_no:gnd_T4_contact_no+2,gnd_T4_contact_no:gnd_T4_contact_no+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,gnd_T4_contact_no:gnd_T4_contact_no+2)=eye(3,3);
body_contact_no=body_contact_no+3;
system_contact_ctr=system_contact_ctr; %system counter does not increase because this point is already included
disp("Frame T4: ");
disp((frame.cm(1:3)+Ai_body*frame.T4)');



%top_hammer

top_hammer.cp_TH_drill_rod_imp=[-0.07;-0.100;0]; % Bit in impact position
top_hammer.cp_TH_drill_rod_ext=[-0.07;-0.156;0]; % Bit in extended position
Gi_body=Gi_mat(top_hammer.cm(4),top_hammer.cm(5),top_hammer.cm(6));
Ai_body=Ai_mat(top_hammer.cm(4),top_hammer.cm(5),top_hammer.cm(6));
body_index=top_hammer.body_index;
u_i_bar_body=top_hammer.cp_TH_drill_rod_imp;

J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;

R_mbd(system_contact_ctr:system_contact_ctr+2,system_contact_ctr:system_contact_ctr+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=eye(3,3);
body_contact_no=body_contact_no+3;
TH_system_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;

disp("Top Hammer: ");
disp((top_hammer.cm(1:3)+Ai_body*u_i_bar_body)');


%piston
piston.cp_piston_tool=[0;-0.2;0]; %contact point between piston and tool
Gi_body=Gi_mat(piston.cm(4),piston.cm(5),piston.cm(6));
Ai_body=Ai_mat(piston.cm(4),piston.cm(5),piston.cm(6));
body_index=piston.body_index;
u_i_bar_body=piston.cp_piston_tool;

J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;

R_mbd(system_contact_ctr:system_contact_ctr+2,system_contact_ctr:system_contact_ctr+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=eye(3,3);
body_contact_no=body_contact_no+3;
piston_contact_no=system_contact_ctr;
system_contact_ctr=system_contact_ctr+3;
disp("Piston: ");
disp((piston.cm(1:3)+Ai_body*u_i_bar_body)');


%drill bit rod
Gi_body=Gi_mat(drill_rod_bit.cm(4),drill_rod_bit.cm(5),drill_rod_bit.cm(6));
Ai_body=Ai_mat(drill_rod_bit.cm(4),drill_rod_bit.cm(5),drill_rod_bit.cm(6));
body_index=drill_rod_bit.body_index;
u_i_bar_body=drill_rod_bit.cp_tool_piston;

J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,TH_system_contact_no:TH_system_contact_no+2)=-eye(3,3);
body_contact_no=body_contact_no+3;

disp("Drill_rod_bit.cp_tool_piston: ");
disp((drill_rod_bit.cm(1:3)+Ai_body*u_i_bar_body)');


u_i_bar_body=drill_rod_bit.cp_tool_top_hammer(:,1);
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
A_mbd(body_contact_no:body_contact_no+2,piston_contact_no:piston_contact_no+2)=-eye(3,3);
body_contact_no=body_contact_no+3;
disp("Drill_rod_bit.cp_tool_top_hammer: ");
disp((drill_rod_bit.cm(1:3)+Ai_body*u_i_bar_body)');

cp=1;
u_i_bar_body=drill_rod_bit.cp_face_button(:,cp);
J_mbd(6*body_index-5:6*body_index-3, body_contact_no:body_contact_no+2)=eye(3,3);
J_mbd(6*body_index-2:6*body_index, body_contact_no:body_contact_no+2)=-Gi_body'*skew_sym(Ai_body*u_i_bar_body)';
Li_mat(body_contact_no:body_contact_no+2,6*body_index-5:6*body_index-3)=eye(3,3);
Li_mat(body_contact_no:body_contact_no+2,6*body_index-2:6*body_index)=-skew_sym(Ai_body*u_i_bar_body)*Gi_body;
R_mbd(system_contact_ctr:system_contact_ctr+2,system_contact_ctr:system_contact_ctr+2)=[Ai_body(:,2),Ai_body(:,1),Ai_body(:,3)]; % Y-axis is normal direction. 
A_mbd(body_contact_no:body_contact_no+2,system_contact_ctr:system_contact_ctr+2)=eye(3,3);
body_contact_no=body_contact_no+3;
system_contact_ctr=system_contact_ctr+3;
disp("Drill_rod_bit.cp_face_button(:,1): ");
disp((drill_rod_bit.cm(1:3)+Ai_body*u_i_bar_body)');

% Construct Bi matrix 

%Cq_sys_dyn is the matrix which contains the elements of Cq matrix for dynamic calculations.
qi=6*Nb-35; % no of independent coordinates

%define position of independent coordinates
qi_pos(1:6)=6*gnd_T1.body_index-5:6*gnd_T1.body_index;
qi_pos(7:12)=6*gnd_T2.body_index-5:6*gnd_T2.body_index;
qi_pos(13:18)=6*gnd_T3.body_index-5:6*gnd_T3.body_index;
qi_pos(19:24)=6*gnd_T4.body_index-5:6*gnd_T4.body_index;
qi_pos(25:30)=6*frame.body_index-5:6*frame.body_index; % x_cm, y_cm, z_cm, phi, theta, psi
qi_pos(31)=6*boom.body_index; % psi
qi_pos(32)=6*rocker.body_index; % psi
qi_pos(33)=6*stick.body_index-1; % theta
qi_pos(34)=6*top_hammer.body_index-4; %y_cm
qi_pos(35)=6*piston.body_index-4; %y_cm
qi_pos(36)=6*rotary_gear.body_index-1; % theta
qi_pos(37)=6*drill_rod_bit.body_index-4; %y_cm


q_trans=row_transpose(q,qi_pos);
Cq_trans=column_transpose(Cq_sys_dyn,qi_pos); 
Cq_i=Cq_trans(:,1:size(qi_pos)(2));
Cq_d=Cq_trans(:,size(qi_pos)(2)+1:size(Cq_trans)(2));
Bi=[eye(13+24,13+24);-inv(Cq_d)*Cq_i];

M_r_tr=row_transpose(M_sys,qi_pos); %rowtranspose
M_rc_tr=column_transpose(M_r_tr,qi_pos); % column transpose of the row transposed matrix

J_mbd_r_tr=row_transpose(J_mbd,qi_pos); %% CHECK JMBD
Li_mat_c_tr=column_transpose(Li_mat,qi_pos);

W_mbd_multiplier=Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd; % W_mbd matrix

W_GCS(3*1-2:3*1,:)=(Li_mat_c_tr(3*5-2:3*5,:)-Li_mat_c_tr(3*1-2:3*1,:))*W_mbd_multiplier; % T1 and gnd
W_GCS(3*2-2:3*2,:)=(Li_mat_c_tr(3*6-2:3*6,:)-Li_mat_c_tr(3*2-2:3*2,:))*W_mbd_multiplier; % T2 and gnd
W_GCS(3*3-2:3*3,:)=(Li_mat_c_tr(3*7-2:3*7,:)-Li_mat_c_tr(3*3-2:3*3,:))*W_mbd_multiplier; % T3 and gnd
W_GCS(3*4-2:3*4,:)=(Li_mat_c_tr(3*8-2:3*8,:)-Li_mat_c_tr(3*4-2:3*4,:))*W_mbd_multiplier; % T4 and gnd
W_GCS(3*5-2:3*5,:)=(Li_mat_c_tr(3*9-2:3*9,:)-Li_mat_c_tr(3*12-2:3*12,:))*W_mbd_multiplier; % TH and drill rod 
W_GCS(3*6-2:3*6,:)=(Li_mat_c_tr(3*10-2:3*10,:)-Li_mat_c_tr(3*11-2:3*11,:))*W_mbd_multiplier; % Piston and drill rod
W_GCS(3*7-2:3*7,:)=Li_mat_c_tr(3*13-2:3*13,:)*W_mbd_multiplier;

W_LCS=R_mbd'*W_GCS*R_mbd;

Cq_sys_dyn_c_tr=column_transpose(Cq_sys_dyn,qi_pos);
D=[J_mbd_r_tr*A_mbd, -Cq_sys_dyn_c_tr'];

G=[(Li_mat_c_tr(3*5-2:3*5,:)-Li_mat_c_tr(3*1-2:3*1,:))
   (Li_mat_c_tr(3*6-2:3*6,:)-Li_mat_c_tr(3*2-2:3*2,:))
   (Li_mat_c_tr(3*7-2:3*7,:)-Li_mat_c_tr(3*3-2:3*3,:))
   (Li_mat_c_tr(3*8-2:3*8,:)-Li_mat_c_tr(3*4-2:3*4,:))
   (Li_mat_c_tr(3*9-2:3*9,:)-Li_mat_c_tr(3*12-2:3*12,:))
   (Li_mat_c_tr(3*10-2:3*10,:)-Li_mat_c_tr(3*11-2:3*11,:))
    Li_mat_c_tr(3*13-2:3*13,:)];
 
R_mbd_gmd =[R_mbd, zeros(21,35)
            zeros(35,21), eye(35,35)];    
    
GMD=Cq_sys_dyn_c_tr*inv(M_rc_tr)*D*R_mbd_gmd;

W_lcs_bar=[W_LCS, zeros(21,35)
           zeros(35,21), zeros(35,35)];
A_eq=GMD;

Vi=zeros(3*7,1); %initial velocities
Vi(3*6-2,1)=-10;
Vi_bar=[Vi;zeros(35,1)];
mu=0*ones(7,1); %coefficient of friction
%mu=[0.3;0.3;0.3;0.3;0;0;0.3];
erest=ones(7,1); %coefficient of restitution
alpha=0.5*ones(7,1); %allocation parameter
%alpha(1:7)=[0.1;0.1;0.1;0.1;0.5;0.5;0.5];
alpha(1:7)=[-0.5;-0.5;-0.5;-0.5;0.5;0.5;0.5];
%alpha(1:7)=[-0.5;-0.5;-0.5;-0.5;0.5;0.5;1];
%alpha(1:7)=[-5;-1;-5;-1;0.5;0.1;1];
%alpha(1:7)=[-5;-1;-5;-1;0.5;0.1;1];
%alpha(1:7)=[-0.5;-0.5;-0.5;-0.5;0.5;0.1;1];
%alpha(1:7)=[-0.2;0.297;0.144;0.2;1.284;1.5;1.5];
%alpha(1:7)=[0;0;0;0;0;0;0]
%P=drem(W_LCS,Vi,mu,erest,alpha); 

P=drem(W_lcs_bar,Vi_bar,A_eq,mu,erest,alpha); 
P(1:21)
disp("Constraint forces and moments");
P(22:56)
disp("Change in velocity in mm / sec");
W_lcs_bar*P*1000
disp("Change in velocity of indeoendent coordinates of the system");
(inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P(1:21))(25:25+12)*1000
%disp("Change in velocity of cm of the linkages only");
del_qi_dot=inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P(1:21)
del_q_dot=Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P(1:21)
lambda=-inv(Cq_sys_dyn_c_tr*Cq_sys_dyn_c_tr')*Cq_sys_dyn_c_tr*(M_rc_tr*del_q_dot-J_mbd_r_tr*A_mbd*R_mbd*P(1:21))

q_sys=(Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P(1:21));

q_dot_trans_mm=arow_transpose(q_sys,qi_pos)(25:24+(8*6));

disp("Change in velocity of frame cm");
q_dot_trans_mm(1:6)
link=frame;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_frame=A'*G*q_dot_trans_mm(4:6);

disp("Change in velocity of boom cm");
q_dot_trans_mm(7:12)
link=boom;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_boom=A'*G*q_dot_trans_mm(10:12);

disp("Change in velocity of rocker cm");
q_dot_trans_mm(13:18)
link=rocker;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_rocker=A'*G*q_dot_trans_mm(16:18);

disp("Change in velocity of stick cm");
q_dot_trans_mm(19:24)
link=stick;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_stick=A'*G*q_dot_trans_mm(22:24);

disp("Change in velocity of top hammer cm");
q_dot_trans_mm(25:30)
link=top_hammer;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_top_hammer=A'*G*q_dot_trans_mm(28:30);

disp("Change in velocity of piston cm");
q_dot_trans_mm(31:36)
link=piston;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_piston=A'*G*q_dot_trans_mm(34:36);

disp("Change in velocity of rotary gear cm");
q_dot_trans_mm(37:42)
link=rotary_gear;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_rotary_gear=A'*G*q_dot_trans_mm(40:42);

disp("Change in velocity of drill rod cm");
q_dot_trans_mm(43:48)
link=drill_rod_bit;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_drill_rod_bit=A'*G*q_dot_trans_mm(46:48);

disp("For comparing with ADAMS:");
q_dot_compare(1:3,1)=q_dot_trans_mm(1:3)*1000;
q_dot_compare(4:6,1)=rad2deg(omega_bar_frame(1:3));
q_dot_compare(7,1)=rad2deg(omega_bar_boom(3));
q_dot_compare(8,1)=rad2deg(omega_bar_rocker(3));
q_dot_compare(9,1)=rad2deg(omega_bar_stick(1));
q_dot_compare(10,1)=q_dot_trans_mm(26)*1000;
q_dot_compare(11,1)=q_dot_trans_mm(32)*1000;
q_dot_compare(12,1)=rad2deg(omega_bar_rotary_gear(2));
q_dot_compare(13,1)=q_dot_trans_mm(44)*1000;
q_dot_compare

%xlswrite('drill results for alpha values.xlsx',q_dot_compare)

%---------COMPARE RESULTS WITH LCP ----------------
function wn = calc_Wn(Bi,J_mbd_r_tr,A_mbd,R_mbd)
  Nc=size(R_mbd)(2)/3;
  W_imp=Bi'*J_mbd_r_tr*A_mbd*R_mbd;
  for i=1:Nc
    wn(:,i)=W_imp(:,3*i-2);
  endfor
endfunction

A=Ai_mat(piston.cm(4),piston.cm(5),piston.cm(6));
q_dot_init_piston=A*[0;Vi(3*6-2,1);0];
q_i_dot_init=zeros(37,1); %24 coordinates for ground and 13 for linkages
q_i_dot_init(24+11)=q_dot_init_piston(2);
e=1*eye(7); %coefficient of restitution
M_LCP=Bi'*M_rc_tr*Bi;
W_N=calc_Wn(Bi,J_mbd_r_tr,A_mbd,R_mbd); % W_N is the matrix which when multiplied with the normal component of impulses results in generalized impulses (Pe) with normal components only. 
gamma_A_N=W_N'*q_i_dot_init;
A=W_N'*inv(M_LCP)*W_N;
B=W_N'*q_i_dot_init+e*gamma_A_N;
[ep,P_LCP_N,retcode]= LCPSolve(A,B,1e-8,1e4);
P_LCP_N
P_LCP=zeros(3*7,1);
P_LCP(1,1)=P_LCP_N(1);
P_LCP(3*2-2,1)=P_LCP_N(2);
P_LCP(3*3-2,1)=P_LCP_N(3);
P_LCP(3*4-2,1)=P_LCP_N(4);
P_LCP(3*5-2,1)=P_LCP_N(5);
P_LCP(3*6-2,1)=P_LCP_N(6);
P_LCP(3*7-2,1)=P_LCP_N(7);

del_qi_dot_LCP=inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P_LCP
del_q_dot_LCP=Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P_LCP
lambda_LCP=-inv(Cq_sys_dyn_c_tr*Cq_sys_dyn_c_tr')*Cq_sys_dyn_c_tr*(M_rc_tr*del_q_dot-J_mbd_r_tr*A_mbd*R_mbd*P_LCP)

q_sys=(Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P_LCP);

q_dot_trans_mm=arow_transpose(q_sys,qi_pos)(25:24+(8*6));

%disp("Change in velocity of frame cm");
q_dot_trans_mm(1:6);
link=frame;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_frame=A'*G*q_dot_trans_mm(4:6);

%disp("Change in velocity of boom cm");
q_dot_trans_mm(7:12);
link=boom;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_boom=A'*G*q_dot_trans_mm(10:12);

%disp("Change in velocity of rocker cm");
q_dot_trans_mm(13:18);
link=rocker;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_rocker=A'*G*q_dot_trans_mm(16:18);

%disp("Change in velocity of stick cm");
q_dot_trans_mm(19:24);
link=stick;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_stick=A'*G*q_dot_trans_mm(22:24);

%disp("Change in velocity of top hammer cm");
q_dot_trans_mm(25:30);
link=top_hammer;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_top_hammer=A'*G*q_dot_trans_mm(28:30);

%disp("Change in velocity of piston cm");
q_dot_trans_mm(31:36);
link=piston;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_piston=A'*G*q_dot_trans_mm(34:36);

%disp("Change in velocity of rotary gear cm");
q_dot_trans_mm(37:42);
link=rotary_gear;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_rotary_gear=A'*G*q_dot_trans_mm(40:42);

%disp("Change in velocity of drill rod cm");
q_dot_trans_mm(43:48);
link=drill_rod_bit;
A=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
G=Gi_mat(link.cm(4),link.cm(5),link.cm(6));
omega_bar_drill_rod_bit=A'*G*q_dot_trans_mm(46:48);

disp("For comparing with ADAMS:");
q_dot_compare(1:3,1)=q_dot_trans_mm(1:3)*1000;
q_dot_compare(4:6,1)=rad2deg(omega_bar_frame(1:3));
q_dot_compare(7,1)=rad2deg(omega_bar_boom(3));
q_dot_compare(8,1)=rad2deg(omega_bar_rocker(3));
q_dot_compare(9,1)=rad2deg(omega_bar_stick(1));
q_dot_compare(10,1)=q_dot_trans_mm(26)*1000;
q_dot_compare(11,1)=q_dot_trans_mm(32)*1000;
q_dot_compare(12,1)=rad2deg(omega_bar_rotary_gear(2));
q_dot_compare(13,1)=q_dot_trans_mm(44)*1000;
q_dot_compare
%{
%beeta=1;
%P_bouncy=(1-beeta-(2*beeta*P'*Vi)/(P'*W_LCS*P))*P;
P=zeros(3*6,1);
P(3*4-2,1)=1000;
P_joint=[-5.125961381
   -4.296696405
    2.756947095
    2.046214892
    2.101650144
   -3.660216214
   -7.276122860
    2.428371547
   -1.879789127
    1.084747383
   -1.984003488
   -7.297300953
    1.933675618
   -5.177936486
   -1.126642412
   54.422894790
   -0.624954517
    4.157240357
    1.554103906
    0.314387960
    0.000748414
   -0.003741793
    0.232397523
    0.857147291
    0.000088487
   -0.286312114
   54.659937172
    5.228240409
    1.391507955
   12.929014283
   -1.860273143
   15.997216852
   -0.486082762
    3.182083519
   -0.000179391];

P_combined=[P;P_joint];
del_q_dot_cp=GMD*P_combined;
%del_q_dot=inv(M_rc_tr)*D*R_mbd_gmd*P_combined;
%del_q_dot= arow_transpose(del_q_dot,qi_pos);
for i=1:3*6
  disp(del_q_dot_cp(i,1));
endfor
%del_q_dot_cp
%{   

del_V_cm=1000*Bi*inv(Bi'*M_rc_tr*Bi)*Bi'*J_mbd_r_tr*A_mbd*R_mbd*P;
del_V_cm= arow_transpose(del_V_cm,qi_pos);
disp("Change in Velocity of cm");
disp("Frame:");
del_V_cm(1:6)
disp("Boom:");
del_V_cm(7:12)
disp("Rocker:");
del_V_cm(13:18)
disp("Stick:");
del_V_cm(19:24)
disp("Top Hammer Housing:");
del_V_cm(25:30)
disp("Pison:");
del_V_cm(31:36)
disp("Rotary Gear:");
del_V_cm(37:42)
disp("Drill Rod Bit:");
del_V_cm(43:48)

Mat_write_csv=[del_V_cm(1:6),del_V_cm(7:12),del_V_cm(13:18),del_V_cm(19:24),del_V_cm(25:30),del_V_cm(31:36),del_V_cm(37:42),del_V_cm(43:48)];
csvwrite('C:\Users\Koushik\MS_IITM\Academics\Research\MBD\Code for 3D mechanisms\impulse.csv',Mat_write_csv);

%% PLOT FIGURE IN 2D plane 
frame_2d.cm=Rot_y(-init_rot)*frame.cm(1:3);
boom_2d.cm=Rot_y(-init_rot)*boom.cm(1:3);
rocker_2d.cm=Rot_y(-init_rot)*rocker.cm(1:3);
stick_2d.cm=Rot_y(-init_rot)*stick.cm(1:3);
top_hammer_2d.cm=Rot_y(-init_rot)*top_hammer.cm(1:3);
drill_rod_bit_2d.cm=Rot_y(-init_rot)*drill_rod_bit.cm(1:3);



Ai_frame=Ai_mat(frame.cm(4),frame.cm(5),frame.cm(6));
frame.T1_GCS=frame.cm(1:3)+Ai_frame*frame.T1;
frame.T2_GCS=frame.cm(1:3)+Ai_frame*frame.T2;
frame.T3_GCS=frame.cm(1:3)+Ai_frame*frame.T3;
frame.T4_GCS=frame.cm(1:3)+Ai_frame*frame.T4;
frame.pt1_GCS=frame.cm(1:3)+Ai_frame*frame.pt1;
frame.pt2_GCS=frame.cm(1:3)+Ai_frame*frame.pt2;

x=[frame.T1_GCS(1),frame.T2_GCS(1),frame.T3_GCS(1),frame.T4_GCS(1),frame.T1_GCS(1)];
y=[frame.T1_GCS(2),frame.T2_GCS(2),frame.T3_GCS(2),frame.T4_GCS(2),frame.T1_GCS(2)];
z=[frame.T1_GCS(3),frame.T2_GCS(3),frame.T3_GCS(3),frame.T4_GCS(3),frame.T1_GCS(3)];

figure;

plot3(x,z,y,"color","green");
x=[frame.T1_GCS(1),frame.cm(1)];
y=[frame.T1_GCS(2),frame.cm(2)];
z=[frame.T1_GCS(3),frame.cm(3)];
hold on;
plot3(x,z,y,"color","green");
x=[frame.T2_GCS(1),frame.cm(1)];
y=[frame.T2_GCS(2),frame.cm(2)];
z=[frame.T2_GCS(3),frame.cm(3)];
plot3(x,z,y,"color","green");
x=[frame.T3_GCS(1),frame.cm(1)];
y=[frame.T3_GCS(2),frame.cm(2)];
z=[frame.T3_GCS(3),frame.cm(3)];
plot3(x,z,y,"color","green");
x=[frame.T4_GCS(1),frame.cm(1)];
y=[frame.T4_GCS(2),frame.cm(2)];
z=[frame.T4_GCS(3),frame.cm(3)];
plot3(x,z,y,"color","green");
x=[frame.pt1_GCS(1),frame.cm(1)];
y=[frame.pt1_GCS(2),frame.cm(2)];
z=[frame.pt1_GCS(3),frame.cm(3)];
plot3(x,z,y,"color","green");
x=[frame.pt2_GCS(1),frame.cm(1)];
y=[frame.pt2_GCS(2),frame.cm(2)];
z=[frame.pt2_GCS(3),frame.cm(3)];
plot3(x,z,y,"color","green");

Ai_boom=Ai_mat(boom.cm(4),boom.cm(5),boom.cm(6));
boom.pt1_GCS=boom.cm(1:3)+Ai_boom*boom.pt1;
boom.pt3_GCS=boom.cm(1:3)+Ai_boom*boom.pt3;
boom.pt5_GCS=boom.cm(1:3)+Ai_boom*boom.pt5;
boom.pt4_GCS=boom.cm(1:3)+Ai_boom*boom.pt4;
x=[boom.pt1_GCS(1),boom.cm(1)];
y=[boom.pt1_GCS(2),boom.cm(2)];
z=[boom.pt1_GCS(3),boom.cm(3)];
plot3(x,z,y,"color","blue");

x=[boom.pt5_GCS(1),boom.cm(1)];
y=[boom.pt5_GCS(2),boom.cm(2)];
z=[boom.pt5_GCS(3),boom.cm(3)];
plot3(x,z,y,"color","blue");


x=[boom.pt3_GCS(1),boom.cm(1)];
y=[boom.pt3_GCS(2),boom.cm(2)];
z=[boom.pt3_GCS(3),boom.cm(3)];
plot3(x,z,y,"color","blue");

x=[boom.pt4_GCS(1),boom.cm(1)];
y=[boom.pt4_GCS(2),boom.cm(2)];
z=[boom.pt4_GCS(3),boom.cm(3)];
plot3(x,z,y,"color","blue");


x=[boom.pt3_GCS(1),boom.pt1_GCS(1)];
y=[boom.pt3_GCS(2),boom.pt1_GCS(2)];
z=[boom.pt3_GCS(3),boom.pt1_GCS(3)];
plot3(x,z,y,"color","blue");

x=[boom.pt3_GCS(1),boom.pt5_GCS(1)];
y=[boom.pt3_GCS(2),boom.pt5_GCS(2)];
z=[boom.pt3_GCS(3),boom.pt5_GCS(3)];
plot3(x,z,y,"color","blue");


x=[boom.pt4_GCS(1),boom.pt1_GCS(1)];
y=[boom.pt4_GCS(2),boom.pt1_GCS(2)];
z=[boom.pt4_GCS(3),boom.pt1_GCS(3)];
plot3(x,z,y,"color","blue");

x=[boom.pt4_GCS(1),boom.pt5_GCS(1)];
y=[boom.pt4_GCS(2),boom.pt5_GCS(2)];
z=[boom.pt4_GCS(3),boom.pt5_GCS(3)];
plot3(x,z,y,"color","blue");

Ai_rocker=Ai_mat(rocker.cm(4),rocker.cm(5),rocker.cm(6));
rocker.pt5_GCS=rocker.cm(1:3)+Ai_rocker*rocker.pt5;
rocker.pt6_GCS=rocker.cm(1:3)+Ai_rocker*rocker.pt6;
rocker.pt7_GCS=rocker.cm(1:3)+Ai_rocker*rocker.pt7;
rocker.pt8_GCS=rocker.cm(1:3)+Ai_rocker*rocker.pt8;

x=[rocker.pt5_GCS(1),rocker.cm(1)];
y=[rocker.pt5_GCS(2),rocker.cm(2)];
z=[rocker.pt5_GCS(3),rocker.cm(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt6_GCS(1),rocker.cm(1)];
y=[rocker.pt6_GCS(2),rocker.cm(2)];
z=[rocker.pt6_GCS(3),rocker.cm(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt7_GCS(1),rocker.cm(1)];
y=[rocker.pt7_GCS(2),rocker.cm(2)];
z=[rocker.pt7_GCS(3),rocker.cm(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt8_GCS(1),rocker.cm(1)];
y=[rocker.pt8_GCS(2),rocker.cm(2)];
z=[rocker.pt8_GCS(3),rocker.cm(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt5_GCS(1),rocker.pt6_GCS(1)];
y=[rocker.pt5_GCS(2),rocker.pt6_GCS(2)];
z=[rocker.pt5_GCS(3),rocker.pt6_GCS(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt6_GCS(1),rocker.pt7_GCS(1)];
y=[rocker.pt6_GCS(2),rocker.pt7_GCS(2)];
z=[rocker.pt6_GCS(3),rocker.pt7_GCS(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt7_GCS(1),rocker.pt8_GCS(1)];
y=[rocker.pt7_GCS(2),rocker.pt8_GCS(2)];
z=[rocker.pt7_GCS(3),rocker.pt8_GCS(3)];
plot3(x,z,y,"color","red");

x=[rocker.pt8_GCS(1),rocker.pt5_GCS(1)];
y=[rocker.pt8_GCS(2),rocker.pt5_GCS(2)];
z=[rocker.pt8_GCS(3),rocker.pt5_GCS(3)];
plot3(x,z,y,"color","red");
