function backazimuth = dsp2bazi_onsite(Z,E,N,interval)

ni     = numel(interval);
covMat = 1./ni*[Z(interval)*Z(interval)', Z(interval)*N(interval)', Z(interval)*E(interval)'; ...
                N(interval)*Z(interval)', N(interval)*N(interval)', N(interval)*E(interval)'; ...
                E(interval)*Z(interval)', E(interval)*N(interval)', E(interval)*E(interval)'];

[eigenVect,eigenVal] = eig(covMat);

% azi = atand(eigenvec[1][0]/eigenvec[2][0]);
azi = atand( eigenVect(2,1)/eigenVect(3,1) );

if eigenVect(3,1)>0; shift = 90;
else                 shift = 270;
end
backazimuth = shift-azi;

if eigenVect(1,1)>0
    if backazimuth<180; backazimuth=backazimuth+180;
    else                backazimuth=backazimuth-180;
    end
end


%    bool calcaz(std::vector<float>& N, std::vector<float>& E, std::vector<float>& Z, d
% ouble &epi_azimuth) {    // Distance, azimuth and coordinates of the epicenter
%        double azi = 0.0;
%        double covmat[3][3];
%        double eigenvec[3][3];
%        double eigenval[3];
%        for (int j=0; j<3; j++) {
%            for (int k=0; k<3; k++) {
%                covmat[j][k] = 0.0;
%            }
%        }
%        for (unsigned int i=0; i<N.size(); i++) {
%            covmat[0][0] += Z[i]*Z[i];
%            covmat[0][1] += Z[i]*N[i];
%            covmat[0][2] += Z[i]*E[i];
%            //covmat[1][0] += N[i]*Z[i];
%            covmat[1][1] += N[i]*N[i];
%            covmat[1][2] += N[i]*E[i];
%            //covmat[2][0] += E[i]*Z[i];
%            //covmat[2][1] += E[i]*N[i];
%            covmat[2][2] += E[i]*E[i];
% 
%            //covmat[0][0] += dis_Z[i]*dis_Z[i];
%            //covmat[0][1] += dis_Z[i]*dis_N[i];
%            //covmat[0][2] += dis_Z[i]*dis_E[i];
%            ////covmat[1][0] += dis_N[i]*dis_Z[i];
%            //covmat[1][1] += dis_N[i]*dis_N[i];
%            //covmat[1][2] += dis_N[i]*dis_E[i];
%            ////covmat[2][0] += dis_E[i]*dis_Z[i];
%            ////covmat[2][1] += dis_E[i]*dis_N[i];
%            //covmat[2][2] += dis_E[i]*dis_E[i];
%        }
%        if (dsyevv3(covmat,eigenvec,eigenval) == -1) {
%            return false;
%        }
%        azi = atand(eigenvec[1][0]/eigenvec[2][0]);
% 
%        int shift = 0;
%        //double epi_azimuth = 0.0;
%        if (eigenvec[2][0] > 0) { // correct for azimuth quadrant
%            shift = 90.0;
%            epi_azimuth = shift - azi;
%        }
%        else {
%            shift = 270.0;
%            epi_azimuth = shift - azi;
%        }
%        if (eigenvec[0][0] > 0) { // correct 180. azimuthal uncertainty using
%            if (epi_azimuth < 180.) {
%                epi_azimuth += 180.0;
%            }
%            else {
%                epi_azimuth -= 180.0;
%            }
% 
%            //WPLib::printLock();
%            //cout << "geocalc: " << eigenvec[0][0] << " " << eigenvec[1][0] << " " <<
% eigenvec[2][0] << " " << azi << " " << epi_azimuth << endl;
%            //cout <<"geocalc: "<<dis_N[0]<<" , "<<dis_E[0]<<endl;
%            //WPLib::printUnlock();
%        }
% 
%        return true;
%    }
% **************************************
