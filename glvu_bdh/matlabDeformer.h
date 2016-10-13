#pragma once

#include "deformer.h"
#include "matlab_utils.h"
#include <thread>
#include <chrono>




const char *solver_names[] = { "CVX", "Direct Mosek" };
int numInterpFrames = 200;

extern bool showATB;
extern int viewport[4];
void display();
void loadP2PConstraints();


struct MatlabDeformer : public Deformer
{
    //std::string meshNameSuffix;
    //std::string matlabNamespace("_BDH");
    //float k4bdh = 0.f;
    //bool useLastFz = true;

    //inline std::string matName(const std::string& s) { return s + matlabNamespace; }

    bool updateMatlabPlot;


    int solver;
    float cage_offset;
    int nVirtualVertex;
    int nFixedSample;
    int nDenseEvaluationSample;
    int nActSetPoolSample;
    float p2p_weight;
    float sigma1_upper_bound;
    float sigma2_lower_bound;
    float k_upper_bound;
    bool solver_output;

    bool binarySearchValidMap;
    float timet;

	unsigned int interpAlgorithm;




    bool initiated;

    Eigen::MatrixXcd C; // Cauchy Coordiantes for vertices

    std::vector<int> p2pIdxs;
    std::vector<float> p2pCoords;

    MyMesh &M;

    MatlabDeformer(MatlabDeformer&) = delete;

    MatlabDeformer(MyMesh &m) :M(m), initiated(false),
        updateMatlabPlot(false), solver(1), cage_offset(2e-2f), 
        nVirtualVertex(1), nFixedSample(1), nActSetPoolSample(2000), 
        p2p_weight(100.f), sigma1_upper_bound(7.f), sigma2_lower_bound(0.35f), k_upper_bound(0.8f),
		solver_output(false), binarySearchValidMap(true), timet(0), 
        interpAlgorithm(1) {

        using deformerptr = MatlabDeformer*;

        TwBar *bar = TwNewBar("MatlabDeformer");

        TwDefine(" MatlabDeformer size='220 380' color='255 0 255' text=dark alpha=128 position='5 380' label='BDH Deformer'"); // change default tweak bar size and color
        //TwAddVarRW(bar, "Mesh NameSuffix", TW_TYPE_STDSTRING, &meshNameSuffix, " ");
        //TwAddVarRW(bar, "k", TW_TYPE_FLOAT, &k4bdh, " min=0 max=1 step=0.1 ");
        //TwAddVarRW(bar, "Update Matlab plot", TW_TYPE_BOOLCPP, &updateMatlabPlot, " ");
        //TwAddVarRW(bar, "Use last fz", TW_TYPE_BOOLCPP, &useLastFz, " ");

        //////////////////////////////////////////////////////////////////////////
        //TwAddSeparator(bar, " Parameters initialization ", " ");

        TwAddVarCB(bar, "cage offset", TW_TYPE_FLOAT,
            [](const void *v, void *d) {  deformerptr(d)->cage_offset = *(const float*)(v); deformerptr(d)->preprocess(); },
            [](void *v, void *d)       { *(float*)(v) = deformerptr(d)->cage_offset; },
            this, " ");

        TwAddVarCB(bar, "numVirtualVertex", TW_TYPE_INT32,
            [](const void *v, void *d) { deformerptr(d)->nVirtualVertex = *(const int*)(v); deformerptr(d)->preprocess(); },
            [](void *v, void *d)       { *(int*)(v) = deformerptr(d)->nVirtualVertex; },
            this, " ");

        TwAddVarCB(bar, "numFixedSample", TW_TYPE_INT32,
            [](const void *v, void *d) { deformerptr(d)->nFixedSample = *(const int*)(v); deformerptr(d)->preprocess(); },
            [](void *v, void *d)       { *(int*)(v) = deformerptr(d)->nFixedSample; },
            this, " ");

        TwAddVarCB(bar, "numDenseEvaluationSample", TW_TYPE_INT32,
            [](const void *v, void *d) { deformerptr(d)->nDenseEvaluationSample = *(const int*)(v); deformerptr(d)->preprocess(); },
            [](void *v, void *d)       { *(int*)(v) = deformerptr(d)->nDenseEvaluationSample; },
            this, " ");

        TwAddVarCB(bar, "numActiveSetPoolSample", TW_TYPE_INT32,
            [](const void *v, void *d) { deformerptr(d)->nActSetPoolSample = *(const int*)(v);  deformerptr(d)->preprocess(); },
            [](void *v, void *d)       { *(int*)(v) = deformerptr(d)->nActSetPoolSample; },
            this, " ");

        //////////////////////////////////////////////////////////////////////////
        TwAddVarRW(bar, "Sigma1", TW_TYPE_FLOAT, &sigma1_upper_bound, " group=Deformer ");
        TwAddVarRW(bar, "Sigma2", TW_TYPE_FLOAT, &sigma2_lower_bound, " group=Deformer ");
        TwAddVarRW(bar, "k", TW_TYPE_FLOAT, &k_upper_bound, " group=Deformer ");

        TwType solverType = TwDefineEnumFromString("Solver", "CVX, Direct Mosek");
        TwAddVarRW(bar, "Solver", solverType, &solver, " group=Deformer help='Least Square is only for conformal map(k=0)' ");

        //TwAddVarRW(bar, "#Iteration", TW_TYPE_INT32, &numIterationPerDeform, " min=1 ");
        TwAddButton(bar, "Reset Mesh", [](void *d){
			matlabEval(" Phi = offsetCage; Psy = Phi * 0; rot_trans=[1; 0];");
            deformerptr(d)->deformResultFromMaltab("PhiPsy"); }, this, " group=Deformer ");


        //////////////////////////////////////////////////////////////////////////
        //TwAddSeparator(bar, " Interpolation ", " ");

        TwAddButton(bar, "Add Keyframe", [](void *){
			matlabEval("PhiPsyKF(:,end+(1:2)) = [Phi Psy];"); }, nullptr, " group=Interpolator ");

        TwAddButton(bar, "Set as Keyframe", [](void *){
			matlabEval("PhiPsyKF(:,ikeyframe*2+(-1:0)) = [Phi Psy];"); }, nullptr, " group=Interpolator ");

        TwAddButton(bar, "Save all Keyframes", [](void *){
			matlabEval("save([datadir 'PhiPsyKF'], 'PhiPsyKF');"); }, nullptr, " group=Interpolator ");

        TwAddVarCB(bar, "View Keyframe", TW_TYPE_INT32,
			[](const void *v, void *d) {
	        scalar2matlab("ikeyframe", *(const int*)(v)); 
			matlabEval("Phi = PhiPsyKF(:, ikeyframe*2-1);  Psy = PhiPsyKF(:, ikeyframe*2);");
            matlabEval("rot_trans=[1; 0];");
			deformerptr(d)->deformResultFromMaltab("PhiPsy");


			matlabEval("P2PCurrentPositions = C(P2PVtxIds, :)*Phi + conj(C(P2PVtxIds, :)*Psy);");
			//loadP2PConstraints();
		},
            [](void *v, void *)       { *(int*)(v) = matlab2scalar("ikeyframe"); },
            this, " min=1 max=1000 step=1 group=Interpolator ");


        TwAddVarCB(bar, "t", TW_TYPE_FLOAT,
            [](const void *v, void *d) {
            deformerptr(d)->timet = *(const float*)(v);
			deformerptr(d)->interpAtTime();
        },
            [](void *v, void *d)       { *(float*)(v) = deformerptr(d)->timet; },
            this, " min=-0.1 max=1.1 step=0.002 keyincr=RIGHT keydecr=LEFT group=Interpolator ");


		TwType InterpAlg = TwDefineEnumFromString("InterpAlgorithm", "BDH/metric, BDH/eta, BDH/nu, SIG13, BDH_GM");
		TwAddVarCB(bar, "InterpAlg", InterpAlg, 
            [](const void *v, void *d) {
			deformerptr(d)->interpAlgorithm = *(unsigned int*)(v);
			deformerptr(d)->interpAtTime();
        },
            [](void *v, void *d)       { *(unsigned int*)(v) = deformerptr(d)->interpAlgorithm; },
			
			this, " group=Interpolator ");


        TwAddVarRW(bar, "numFrame", TW_TYPE_INT32, &numInterpFrames, " group=Interpolator ");


        TwAddButton(bar, "Generate sequence", [](void *d){
            auto * pthis = deformerptr(d);

            vector2matlab("interpAnchID", pthis->M.auxVtxIdxs);
            matlabEval("interpAnchID= interpAnchID+1;"); 
            matlabEval("bdhInterpPrep");

            auto &M = pthis->M;
            auto constrainVertices0 = M.constrainVertices;
            auto constrainVerticesRef0 = M.constrainVerticesRef;
            auto auxVtxIdxs0 = M.auxVtxIdxs;
            M.auxVtxIdxs.clear();
            M.constrainVertices.clear();
            M.constrainVerticesRef.clear();


            bool showATB0 = showATB;
            float edgeWidth0 = M.edgeWidth;
            //float meshScale0 = M.mMeshScale;
            //auto translate = M.mTranslate;

           
            {

                int vw = viewport[2]/8*8, vh = viewport[3];
                std::vector<unsigned char> pixels(vw*vh * 3);


                for (int i = 0; i < numInterpFrames; i++) {
                    pthis->timet = i / (numInterpFrames - 1.f);
                    pthis->interpAtTime();
                    //std::string imgfile = datadir + "morph_" + std::to_string(i) + ".jpg";

                    display();

                    glReadPixels(viewport[0], viewport[1], vw, vh, GL_BGR, GL_UNSIGNED_BYTE, &pixels[0]);
                }
            }


            M.constrainVertices = constrainVertices0;
            M.constrainVerticesRef = constrainVerticesRef0;
            M.auxVtxIdxs = auxVtxIdxs0;
        }, this, " key=m group=Interpolator ");

       preprocess();
       initiated = true;
    }

    ~MatlabDeformer(){
        TwBar *bar = TwGetBarByName("MatlabDeformer");
        if (bar)    TwDeleteBar(bar); 
    }

    virtual std::string name(){ return "P2PHarmonic"; }

    const char* currentInterpAlgorithm() const
    {
        const char* algNames[] = { "BDHI/metric", "BDHI/eta",  "BDHI/nu",   "Chen13", "ARAP", "ARAP_LG", "FFMP", "BDH_GM", "Unknown" }; // sync with TwType InterpAlg
        return algNames[std::min<int>(std::size(algNames) - 1, interpAlgorithm)];
    }

    virtual void preprocess() 
    {
        if (initiated){
            scalar2matlab("numVirtualVertices", nVirtualVertex);
            scalar2matlab("numFixedSamples", nFixedSample);

            scalar2matlab("cage_offset", cage_offset);
            scalar2matlab("numDenseEvaluationSamples", nDenseEvaluationSample);
            scalar2matlab("numActiveSetPoolSamples", nActSetPoolSample);
        }

        matlabEval("p2p_harmonic_prep;");

        C = matlab2eigenComplex("C");

        nVirtualVertex = (int)matlab2scalar("numVirtualVertices");
        nFixedSample = (int)matlab2scalar("numFixedSamples");

        if (!initiated){
            cage_offset = (float)matlab2scalar("cage_offset");
            nDenseEvaluationSample = (int)matlab2scalar("numDenseEvaluationSamples");
            nActSetPoolSample = (int)matlab2scalar("numActiveSetPoolSamples");

            if ((bool)matlab2scalar("hasDataLoaded")){
                sigma1_upper_bound = (float)matlab2scalar("sigma1_upper_bound");
                sigma2_lower_bound = (float)matlab2scalar("sigma2_lower_bound");
                k_upper_bound = (float)matlab2scalar("k_upper_bound");
            }

            //solver_output = !(bool)matlab2scalar("no_output");
            //p2p_weight = (float)matlab2scalar("p2p_weight");
            //binarySearchValidMap = (bool)matlab2scalar("binarySearchValidMap");
        }

        deformResultFromMaltab("PhiPsy");
    }


	void interpAtTime() {
        float t = timet;

		const char *bdhimethods[] = { "metric", "eta", "nu" };
        
        // We pass in interpAlgorithm to specify which algorithm to use
        // in fBdhInterpX
		switch (interpAlgorithm) {
		case 0:			// BDH / metric
		case 1:			// BDH / eta
		case 2:			// BDH / nu

			matlabEval(std::string("bdhiMethod='")+bdhimethods[interpAlgorithm]+"';");

			// TODO: clear rot_trans mess, can be merged into Phi Psy
			//matlabEval("[Phi, Psy, rot_trans] = fBdhInterp(" + std::to_string(t) + ");");
            matlabEval("XBDHI = fBdhInterpX(" + std::to_string(t) + ");");
			deformResultFromMaltab("XBDHI");
            if (M.vizVtxData) {
                matlabEval("k = fBdhInterpkX(" + std::to_string(t) + ");");
                auto k = matlab2vector<float>("single(k)", true);
                M.setVertexDataViz(k.data());
            }
			break;
		case 3:			// SIG13
			matlabEval("XSIG13 = fSIG13Interp(" + std::to_string(t) + ");");
			deformResultFromMaltab("XSIG13");
			break;
        case 4:         // BDH_GM, BDHI for general map other than harmonic map
			matlabEval("XBDHGM = fGBDHInterp(" + std::to_string(t) + ");");
			deformResultFromMaltab("XBDHGM");
			break;
		}
	}



    virtual void updateP2PConstraints(int) 
    {
        using namespace Eigen;
        const size_t nConstrain = M.constrainVertices.size();
        eigen2matlab("P2PVtxIds", (Map<VectorXi>(M.getConstrainVertexIds().data(), nConstrain) + VectorXi::Ones(nConstrain)).cast<double>());
        matlabEval("CauchyCoordinatesAtP2Phandles = C(P2PVtxIds,:);");

        MatrixX2d p2pdst = Map<Matrix<float, Dynamic, 2, RowMajor> >(M.getConstrainVertexCoords().data(), nConstrain, 2).cast<double>();
        eigen2matlabComplex("P2PCurrentPositions", p2pdst.col(0), p2pdst.col(1));
#if 0
        eigen2matlab(matName("handleIds"), (Map<VectorXi>(M.getConstrainVertexIds().data(), nConstrain) + VectorXi::Ones(nConstrain)).cast<double>());
        eigen2matlab(matName("handlePos"), Map<MatrixXf>(M.getConstrainVertexCoords().data(), nConstrain, 2).cast<double>());

        bool updatePlotInMatlab = true;
        if (updatePlotInMatlab){
            matlabEval("fUpdatePlotVertices(uihandle.hconstrains," + matName("handlePos") + ");\n"
                + "fUpdatePlotVertices(uihandle.hconstraint,uihandle.hmt.Vertices(" + matName("handleIds") + ",:));\n"
                + "fUpdatePlotVertices(hcs,fC2R(z(" + matName("handleIds") + ")));");
        }
#endif
    }


    void deformResultFromMaltab(std::string resVarName)
    {
        using namespace Eigen;
		if ( !resVarName.compare("PhiPsy") ) { // do all interpolation computation in matlab, for better performance with # virtual vertex > 1
			MatrixXcd Phi = matlab2eigenComplex("Phi");
			MatrixXcd Psy = matlab2eigenComplex("Psy");

			if (Phi.rows() == 0 || Psy.rows() == 0 || C.rows() == 0) return;

			//using Vec = Eigen::Map < Eigen::VectorXcd > ;
			//Eigen::VectorXcd x = C*Vec(Phi.data(), Phi.rows()) + (C*Vec(Psy.data(), Psy.rows())).conjugate();
			Eigen::VectorXcd x = C*Phi + (C*Psy).conjugate();

			if (getMatEngine().hasVar("rot_trans")) {
				// for interpolation
				Vector2cd rot_trans = matlab2eigenComplex("rot_trans");
				x = x.array()*rot_trans(0) + rot_trans(1);
			}

			if (x.rows() == 0) return;

			Matrix<float, Dynamic, 2, RowMajor> xr(x.rows(), 2);
			xr.col(0) = x.real().cast<float>();
			xr.col(1) = x.imag().cast<float>();
			M.upload(xr, Eigen::MatrixXi(), nullptr);
		}
		else {
			MatrixXcd x = matlab2eigenComplex(resVarName);
			if (x.rows() == 0) return;

			Matrix<float, Dynamic, 2, RowMajor> xr(x.rows(), 2);
			xr.col(0) = x.real().cast<float>();
			xr.col(1) = x.imag().cast<float>();
			M.upload(xr, Eigen::MatrixXi(), nullptr);
		}
    }

    virtual void deform()
    {
        //matlabEval(std::string("solver_type='") + solver_names[solver] + "';"  );
        string2matlab("solver_type", solver_names[solver]);
        scalar2matlab("no_output", !solver_output);
        scalar2matlab("p2p_weight", p2p_weight);
        scalar2matlab("sigma1_upper_bound", sigma1_upper_bound);
        scalar2matlab("sigma2_lower_bound", sigma2_lower_bound);
        scalar2matlab("k_upper_bound", k_upper_bound);
        scalar2matlab("binarySearchValidMap", binarySearchValidMap);

        matlabEval("p2p_harmonic;");
        matlabEval("clear rot_trans;");

        deformResultFromMaltab("PhiPsy");

#if 0
        bool updatePlotInMatlab = false;
        scalar2matlab(matName("k"), k4bdh);

        std::stringstream ss;
        if (useLastFz){
            ss << "if exist('" << matName("fz") << "');" << "data.fz=" << matName("fz") << "; end;\n";
        }
        else{
            ss << "if isfield(data,'fz'); data = rmfield(data, 'fz'); end;\n";
        }

        ss << "data.solver=" << solver << ";\n";
        ss << "data.nIter=" << numIterationPerDeform << ";\n";
        matlabEval(ss.str());
        ss.str("");

        ss << "[" << matName("zf") << "," << matName("fz") << "] = BDHarmonicDeform2(data," << matName("k") << ", " << matName("handleIds") << ", " << matName("handlePos")
            << (updatePlotInMatlab ? ",uihandle);" : "); ");
        matlabEval(ss.str());

        Eigen::MatrixXd x;
        matlab2eigen("fC2R(" + matName("zf") + ")", x, true);

        if (x.count() == 0) return;

        M.upload(x.cast<float>().eval(), Eigen::MatrixXi(), nullptr);
#endif
    }

    virtual void resetDeform() {
		matlabEval("Phi = offsetCage; Psy = Phi * 0; clear rot_trans;");
        deformResultFromMaltab("PhiPsy"); 
    }
    virtual void getResult() {}
    virtual void saveData()   { matlabEval("p2p_harmonic_savedata;"); }
};


