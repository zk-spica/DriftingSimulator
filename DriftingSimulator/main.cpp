
#include <cstdio>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>

#include "ogldev_engine_common.h"
#include "ogldev_util.h"
#include "ogldev_app.h"
#include "ogldev_pipeline.h"
#include "ogldev_camera.h"
#include "ogldev_glut_backend.h"
#include "mesh.h"
#include "skybox.h"
#include "ogldev_atb.h"
#include "ogldev_shadow_map_fbo.h"
#include "lighting_technique.h"
#include "shadow_map_technique.h"

#define WINDOW_WIDTH  1024
#define WINDOW_HEIGHT 768

#define abs(x) ((x)>0?(x):(-(x)))
#define sqr(x) ((x)*(x))

#define debug if(false)


static const double wheelAcceleration = 0.006, maxWheelVelocity = 1, wheelDeceleration = 0.01, maxWheelReverseVelocity = 0.3, 
                    turningVelocity = 0.018, velocityEps = 1e-6, bounceFactor = 2; 
static double carMass = 500, carAngularVelocity = 0, wheelVelocity = 0, wheelAngle = 0, combinedAngularAcceleration, wheelRoll,
              maxWheelAngle = 0.37, frictionFactor = 0.01;
static Vector2f carVelocity, combinedForce;
static const float groundHeight = -0.329f;
double carLength = 11, carWidth = 4.5;
float viewYAngle;
float carScale_s = 0.0008f, carScale = 0.0007f;
static double cameraRotateH = 0, cameraRotateV = 0, cameraDis = 15;
const char *cityModelPath = "../Content/venice/venice_0.obj";

bool displaySkyBox = true;
bool displayCar = true;

static float FieldDepth = 10.0f;

Vector3f vertex[1000000];
int nVertex = 0;
float eps3 = 0.001f;

struct vecExample{
	Vector3f vecLoc;
	int isPureColor;

	vecExample(){
		isPureColor = 1;
	}
};

class Line
{
public:
	float x1, y1, x2, y2;
	Line* next;

	Line(float _x1, float _y1, float _x2, float _y2)
	{
		x1 = _x1, y1 = _y1;
		x2 = _x2, y2 = _y2;
		next = NULL;
	}
};

double det(double x1, double y1, double x2, double y2)
{
	return x1*y2 - x2*y1;
}

double dot(double x1, double y1, double x2, double y2)
{
	return x1*x2 + y1*y2;
}


double det(const Vector2f &a, const Vector2f &b)
{
	return a.x*b.y - a.y*b.x;
}

class Face
{
public:
	int *v, nv;
	double A, B, C, D;
	int projectionDir;  //0:xy 1:yz 2:xz
	Face* next;

	Face(int p[], int np)
	{
		v = (int*)malloc((np+1)*sizeof(int));
		memcpy(v, p, (np+1)*sizeof(int));
		nv = np;
		A = det(vertex[v[1]].y-vertex[v[0]].y, vertex[v[1]].z-vertex[v[0]].z, vertex[v[2]].y-vertex[v[0]].y, vertex[v[2]].z-vertex[v[0]].z);
		B = det(vertex[v[1]].z-vertex[v[0]].z, vertex[v[1]].x-vertex[v[0]].x, vertex[v[2]].z-vertex[v[0]].z, vertex[v[2]].x-vertex[v[0]].x);
		C = det(vertex[v[1]].x-vertex[v[0]].x, vertex[v[1]].y-vertex[v[0]].y, vertex[v[2]].x-vertex[v[0]].x, vertex[v[2]].y-vertex[v[0]].y);
		D = -A*vertex[v[0]].x-B*vertex[v[0]].y-C*vertex[v[0]].z;
		next = NULL;

		double s = 0;
		int i;
		for (i=0; i<nv; ++i) s += det(vertex[v[i]].x-vertex[v[0]].x, vertex[v[i]].y-vertex[v[0]].y, vertex[v[i+1]].x-vertex[v[0]].x, vertex[v[i+1]].y-vertex[v[0]].y);
		if (abs(s) > eps3) projectionDir = 0;
		else
		{
			s = 0;
			for (i=0; i<nv; ++i) s += det(vertex[v[i]].y-vertex[v[0]].y, vertex[v[i]].z-vertex[v[0]].z, vertex[v[i+1]].y-vertex[v[0]].y, vertex[v[i+1]].z-vertex[v[0]].z);
			if (abs(s) > eps3) projectionDir = 1; else projectionDir = 2;
		}
	}
};

const float cityMaxX = 3466, cityMaxZ = 1986, cityMinX = 2411, cityMinZ = 1235;
//(X: 2412-3466, Z:1236-1986)
const float chunkWidth = 50;
#define MAXCHUNK 100
Line *chunkLines[MAXCHUNK][MAXCHUNK];
Face *chunkFaces[MAXCHUNK][MAXCHUNK];
int vis[100][100], tag;

void warning()
{
	printf("warning\n");
	for (;;);
}

struct tankModel{
	Vector3f center;
	float yAngle;

	tankModel(){
		center = Vector3f(0.0f , 0.0f, 0.0f);
		yAngle = 0.0f;
	}

};

extern int keyStatus[300];

class Tutorial22 : public ICallbacks
{
public:
    Tutorial22()
    {
        m_pGameCamera = NULL;
        m_pEffect = NULL;
		m_pSkyBox = NULL;

        m_scale = 0.0f;
        m_directionalLight.Color = Vector3f(1.0f, 1.0f, 1.0f);
        m_directionalLight.AmbientIntensity = 0.3f;
        m_directionalLight.DiffuseIntensity = 0.8f;
        m_directionalLight.Direction = Vector3f(-1.0f, -1.0, -1.0);

        m_persProjInfo.FOV = 60.0f;
        m_persProjInfo.Height = WINDOW_HEIGHT;
        m_persProjInfo.Width = WINDOW_WIDTH;
        m_persProjInfo.zNear = 1.0f;
        m_persProjInfo.zFar = 1000.0f;  

		/***********SpotLight*****/
		m_spotlight.DiffuseIntensity = 0.9f;
		m_spotlight.Color = Vector3f(1.0f, 1.0f, 1.0f);
		m_spotlight.Attenuation.Linear = 0.1f;
		m_spotlight.Cutoff = 10.0f;
		/*************/

		/****** Shadow Projection*************/
		m_shadowOrthoProjInfo.l = -80.0f;
		m_shadowOrthoProjInfo.r = 80.0f;
		m_shadowOrthoProjInfo.t = 80.0f;
		m_shadowOrthoProjInfo.b = -80.0f;
		m_shadowOrthoProjInfo.n = -70.0f;
		m_shadowOrthoProjInfo.f = 50.0f;
		/****************************************/
    }

    ~Tutorial22()
    {
        delete m_pEffect;
        delete m_pGameCamera;
        delete m_pMesh;
		delete m_pSkyBox;
    }

    bool Init()
    {
		/* ShadowMap Program and Shadow Frame Buffer */
		if (!m_shadowMapFBO.Init(WINDOW_WIDTH, WINDOW_HEIGHT)) {
			return false;
		}
		if (!m_ShadowMapEffect.Init()) {
			printf("Error initializing the shadow map technique\n");
			return false;
		}
		/***************************************************/

        Vector3f Pos(3007.0f, 10.0f, 1609.0f);
        Vector3f Target(3.0f, -1.0f, 3.0f);
        Vector3f Up(0.0, 1.0f, 0.0f);

        m_pGameCamera = new Camera(WINDOW_WIDTH, WINDOW_HEIGHT, Pos, Target, Up);

        m_pEffect = new LightingTechnique();

        if (!m_pEffect->Init()) {
            printf("Error initializing the lighting technique\n");
            return false;
        }

        m_pEffect->Enable();

		/* Initialize Lighting Program */
		m_pEffect->SetColorTextureUnit(COLOR_TEXTURE_UNIT_INDEX);
		m_pEffect->SetShadowMapTextureUnit(SHADOW_TEXTURE_UNIT_INDEX);
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);

        m_pMesh = new Mesh();
		car_pMesh = new Mesh();
		rearWheel_pMesh = new Mesh();
		frontLeftWheel_pMesh = new Mesh();
		frontRightWheel_pMesh = new Mesh();
		frontLeftBrake_pMesh = new Mesh();
		frontRightBrake_pMesh = new Mesh();


		printf("load car\n");
		car_pMesh->LoadMesh("../Content/BMW_M3_GTR/BMW_M3_GTR.obj");
		rearWheel_pMesh->LoadMesh("../Content/BMW_M3_GTR/rear_wheel.obj");
		frontLeftWheel_pMesh->LoadMesh("../Content/BMW_M3_GTR/front_left_wheel.obj");
		frontRightWheel_pMesh->LoadMesh("../Content/BMW_M3_GTR/front_right_wheel.obj");
		frontLeftBrake_pMesh->LoadMesh("../Content/BMW_M3_GTR/front_left_brake.obj");;
		frontRightBrake_pMesh->LoadMesh("../Content/BMW_M3_GTR/front_right_brake.obj");;

		car.center = Vector3f(3007.0f, -0.00f, 1609.0f);
		//car.center = Vector3f(3312.0f, -0.088314f, 1533.0f);
		

		m_pSkyBox = new SkyBox(m_pGameCamera, m_persProjInfo);


		if (!m_pSkyBox->Init(".",
			"../Content/sp3right.jpg",
			"../Content/sp3left.jpg",
			"../Content/sp3top.jpg",
			"../Content/sp3bot.jpg",
			"../Content/sp3front.jpg",
			"../Content/sp3back.jpg")) {
			return false;
		}

		initializeMapBorder();


		printf("load city\n");
        return m_pMesh->LoadMesh(cityModelPath);
    }

	void getChunk(float x, float y, int &X, int &Y)
	{
		X = floor((x-cityMinX) / chunkWidth);
		Y = floor((y-cityMinZ) / chunkWidth);
		if (X<0 || X>=MAXCHUNK || Y<0 || Y>=MAXCHUNK) X = Y = MAXCHUNK-1;
	}

	double determinant(double v1, double v2, double v3, double v4)
	{
		return v1*v4-v2*v3;
	}

	bool intersect_LineLine(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
	{
		double delta = determinant(x2-x1, x3-x4, y2-y1, y3-y4);
		if ( delta<=(1e-5) && delta>=-(1e-5))
		{
			return false;
		}
		double namenda = determinant(x3-x1, x3-x4, y3-y1, y3-y4) / delta;
		if ( namenda>1+1e-3 || namenda<-1e-3 )
		{
			return false;
		}
		double miu = determinant(x2-x1, x3-x1, y2-y1, y3-y1) / delta;
		if ( miu>1+1e-3 || miu<-1e-3 )
		{
			return false;
		}
		return true;
	}

	bool inChunk(int i, int j, float x, float y)
	{
		return cityMinX+i*chunkWidth-eps3 <= x && x<= cityMinX+(i+1)*chunkWidth+eps3 && cityMinZ+j*chunkWidth-eps3 <= y && y<= cityMinZ+(j+1)*chunkWidth+eps3;
	}

	bool intersect_ChunkLine(int i, int j, float x1, float y1, float x2, float y2)
	{
		if (intersect_LineLine(cityMinX+i*chunkWidth, cityMinZ+j*chunkWidth, cityMinX+(i+1)*chunkWidth, cityMinZ+j*chunkWidth, x1, y1, x2, y2)) return true;
		if (intersect_LineLine(cityMinX+i*chunkWidth, cityMinZ+j*chunkWidth, cityMinX+i*chunkWidth, cityMinZ+(j+1)*chunkWidth, x1, y1, x2, y2)) return true;
		if (intersect_LineLine(cityMinX+(i+1)*chunkWidth, cityMinZ+j*chunkWidth, cityMinX+(i+1)*chunkWidth, cityMinZ+(j+1)*chunkWidth, x1, y1, x2, y2)) return true;
		if (intersect_LineLine(cityMinX+i*chunkWidth, cityMinZ+(j+1)*chunkWidth, cityMinX+(i+1)*chunkWidth, cityMinZ+(j+1)*chunkWidth, x1, y1, x2, y2)) return true;
		if (inChunk(i,j,x1,y1) || inChunk(i,j,x2,y2)) return true;
	}

	void addBarrier(float x1, float y1, float x2, float y2)
	{
//		printf("addbarrier %f %f %f %f\n", x1, y1, x2, y2);

		int X1, Y1, X2, Y2, i, j;
		getChunk(x1, y1, X1, Y1);
		getChunk(x2, y2, X2, Y2);
		if (X1 > X2) swap(X1, X2);
		if (Y1 > Y2) swap(Y1, Y2);
		for (i=X1-1; i<=X2+1; ++i)
		for (j=Y1-1; j<=Y2+1; ++j)
		if (intersect_ChunkLine(i, j, x1, y1, x2, y2))
		{
			Line *newLine = new Line(x1,y1,x2,y2);
			newLine->next = chunkLines[i][j];
			chunkLines[i][j] = newLine;
		}
	}

	void checkBuildingEdge(int v1, int v2)
	{
//		printf("checkBuildingEdge %d %d\n", v1, v2);

		if (abs(vertex[v1].y - vertex[v2].y) > 1) return;
		if (vertex[v1].y < groundHeight || vertex[v2].y < groundHeight) return;
		if (vertex[v1].y > groundHeight + carWidth || vertex[v2].y > groundHeight + carWidth) return;
		addBarrier(vertex[v1].x, vertex[v1].z, vertex[v2].x, vertex[v2].z);
	}

	bool inersectWithBarrier(float x1, float y1, float x2, float y2)
	{
//		printf("intersectWidhtBarrier %f %f %f %f\n", x1, y1, x2, y2);

		int X1, Y1, X2, Y2, i, j;
		Line *p;
		getChunk(x1, y1, X1, Y1);
		getChunk(x2, y2, X2, Y2);
		if (X1 > X2) swap(X1, X2);
		if (Y1 > Y2) swap(Y1, Y2);
		for (i=X1; i<=X2; ++i)
		for (j=Y1; j<=Y2; ++j)
		{
			for (p=chunkLines[i][j]; p!=NULL; p=p->next)
			{
				if (intersect_LineLine(p->x1, p->y1, p->x2, p->y2, x1, y1, x2, y2))
				{
					printf("%f %f %f %f intersect %f %f %f %f\n", x1, y1, x2, y2, p->x1, p->y1, p->x2, p->y2);
					return true;
				}
			}
		}
		return false;
	}

	void checkGroundEdge(int v1, int v2)
	{
		if (vertex[v1] < vertex[v2]) swap(v1, v2);
	}

	void addFace(int i, int j, Face* face)
	{
		face->next = chunkFaces[i][j];
		chunkFaces[i][j] = face;
	}

	bool intersect_LineFace(Vector3f a, Vector3f b, Face *f)
	{
		double s1, s2;
		s1 = a.x*f->A + a.y*f->B + a.z*f->C + f->D;
		s2 = b.x*f->A + b.y*f->B + b.z*f->C + f->D;
		return s1*s2 < 0 || abs(s1) < eps3 || abs(s2) < eps3;
	}

	double module(double x, double y)
	{
		return sqrt(x*x+y*y);
	}

	double angle(double x1, double y1, double x2, double y2)
	{
		if (module(x1,y1)<1e-2 || module(x2,y2)<1e-2) return 0;

		double d1 = det(x1, y1, x2, y2);
		double d2 = dot(x1, y1, x2, y2);
		double theta = asin(d1/module(x1,y1)/module(x2,y2));
		if (d2 < 0)
		{
			if (d1 > 0) theta = M_PI-theta; else theta = -M_PI-theta;
		}
		return theta;
	}

	bool insideFace(Face *f, Vector3f dot)
	{
		double s = 0;
		int i;
		if (f->projectionDir == 0)
		{
			for (i=0; i<f->nv; ++i) s += angle(vertex[f->v[i]].x-dot.x, vertex[f->v[i]].y-dot.y, vertex[f->v[i+1]].x-dot.x, vertex[f->v[i+1]].y-dot.y);
		}
		else
		if (f->projectionDir == 1)
		{
			for (i=0; i<f->nv; ++i) s += angle(vertex[f->v[i]].z-dot.z, vertex[f->v[i]].y-dot.y, vertex[f->v[i+1]].z-dot.z, vertex[f->v[i+1]].y-dot.y);
		}
		else
		if (f->projectionDir == 2)
		{
			for (i=0; i<f->nv; ++i) s += angle(vertex[f->v[i]].x-dot.x, vertex[f->v[i]].z-dot.z, vertex[f->v[i+1]].x-dot.x, vertex[f->v[i+1]].z-dot.z);
		}
		return abs(s) > 1;
	}

	void checkIntersectionWithFace(Vector3f camTar, Vector3f &camPos)
	{
		int X1, Y1, X2, Y2, i, j;
		Face *p;
		double minMiu = 1e10;
		Vector3f dir = camPos-camTar;
		Vector3f Right = Vector3f(0.0f, 1.0f, 0.0f).Cross(camTar-camPos);
        Right.Normalize();
		Right *= WINDOW_WIDTH/853.3;
		Vector3f camPosL, camPosR, camPosD;
		camPosD = camPos-Vector3f(0.0f,WINDOW_HEIGHT/501.8,0.0f);
		camPosL = camPos-Right;
		camPosR = camPos+Right;
		
		Vector3f dirL = camPosL-camTar, dirR = camPosR-camTar, dirD = camPosD-camTar;
		getChunk(camTar.x, camTar.z, X1, Y1);
		getChunk(camPos.x, camPos.z, X2, Y2);
		if (X1 > X2) swap(X1, X2);
		if (Y1 > Y2) swap(Y1, Y2);
		for (i=X1; i<=X2; ++i)
		for (j=Y1; j<=Y2; ++j)
		{
			for (p=chunkFaces[i][j]; p!=NULL; p=p->next)
			{
				if (intersect_LineFace(camTar, camPos, p))
				{
					double miu = (-p->D - p->A*camTar.x - p->B*camTar.y - p->C*camTar.z) / (p->A*dir.x + p->B*dir.y + p->C*dir.z);
					Vector3f dot = camTar+dir*miu;
					if (insideFace(p, dot))
					{
						if (miu < minMiu) minMiu = miu;
					}
				}

				if (intersect_LineFace(camTar, camPosL, p))
				{
					double miu = (-p->D - p->A*camTar.x - p->B*camTar.y - p->C*camTar.z) / (p->A*dirL.x + p->B*dirL.y + p->C*dirL.z);
					Vector3f dot = camTar+dirL*miu;
					if (insideFace(p, dot))
					if (miu < minMiu) minMiu = miu;
				}

				if (intersect_LineFace(camTar, camPosR, p))
				{
					double miu = (-p->D - p->A*camTar.x - p->B*camTar.y - p->C*camTar.z) / (p->A*dirR.x + p->B*dirR.y + p->C*dirR.z);
					Vector3f dot = camTar+dirR*miu;
					if (insideFace(p, dot))
					if (miu < minMiu) minMiu = miu;
				}

				if (intersect_LineFace(camTar, camPosD, p))
				{
					double miu = (-p->D - p->A*camTar.x - p->B*camTar.y - p->C*camTar.z) / (p->A*dirD.x + p->B*dirD.y + p->C*dirD.z);
					Vector3f dot = camTar+dirD*miu;
					if (insideFace(p, dot))
					if (miu < minMiu) minMiu = miu;
				}
			}
		}

		if (minMiu < 1e10-1)
		{
			camPos = camTar+dir*(minMiu-0.001);
		}
//		printf("u=%d\n", u);
	}

	void initializeMapBorder()
	{
		memset(chunkLines, 0, sizeof(chunkLines));
		memset(chunkFaces, 0, sizeof(chunkFaces));

		printf("cityModelPath=%s\n", cityModelPath);

		FILE *fp = fopen(cityModelPath, "r");
		int ch, nPolygonVertex, i;
		int polygonVertex[100];
//		mapGroundEdge.clear();
//		cityMaxX = cityMaxZ = -1e10;
//		cityMinX = cityMinZ = 1e10;
		while (1)
		{
			for (ch = getc(fp); ch <= 32 && ch != -1; ch = getc(fp));
			if (ch == -1) break;

//			printf("ch=%d\n", ch);

			if (ch == 'v')
			{
				ch = getc(fp);
				if (ch == 'n' || ch == 't')
				{
					for (ch = getc(fp); ch !='\n' && ch != '\r'; ch = getc(fp));
					continue;
				}
				++nVertex;
				if (nVertex > 990000) warning();
				fscanf(fp, "%f%f%f", &vertex[nVertex].x, &vertex[nVertex].y, &vertex[nVertex].z);
//				if (vertex[nVertex].x > cityMaxX) cityMaxX = vertex[nVertex].x;
//				if (vertex[nVertex].x < cityMinX) cityMinX = vertex[nVertex].x;
//				if (vertex[nVertex].z > cityMaxZ) cityMaxZ = vertex[nVertex].z;
//				if (vertex[nVertex].z < cityMinZ) cityMinZ = vertex[nVertex].z;
			}
			else
			if (ch == 'f')
			{
				int v1,v2,v3;
				nPolygonVertex = 0;
				while (1)
				{
					fscanf(fp, "%d/%d/%d", &v1, &v2, &v3);
					polygonVertex[++nPolygonVertex] = v1;
					ch = getc(fp);
					if (ch == '\n' || ch == '\r') break;
				}

//				printf("npolygonxvertex=%d\n", nPolygonVertex);

				polygonVertex[0] = polygonVertex[nPolygonVertex];
				bool isGround = true, aboveGround = false;

//				printf("find a surface: "); for (i=1; i<=nPolygonVertex; ++i) printf("%d ", polygonVertex[i]); printf("\n");

				for (i=1; i<=nPolygonVertex; ++i)
				{
					if (abs(vertex[polygonVertex[i]].y - groundHeight) > eps3) isGround = false;
					if (vertex[polygonVertex[i]].y > groundHeight+eps3) aboveGround = true;
				}
				Face *newFace = new Face(polygonVertex, nPolygonVertex);
				++tag;
				for (i = 0; i<nPolygonVertex; ++i)
				{
					if (!isGround && aboveGround) checkBuildingEdge(polygonVertex[i], polygonVertex[i+1]);
					if (isGround) checkGroundEdge(polygonVertex[i], polygonVertex[i+1]);
					int X1, Y1, X2, Y2;

					getChunk(vertex[polygonVertex[i]].x, vertex[polygonVertex[i]].z, X1, Y1);
					getChunk(vertex[polygonVertex[i+1]].x, vertex[polygonVertex[i+1]].z, X2, Y2);
					if (X1 > X2) swap(X1, X2);
					if (Y1 > Y2) swap(Y1, Y2);
					for (int i1=X1; i1<=X2; ++i1)
					for (int j1=Y1; j1<=Y2; ++j1)
					if (vis[i1][j1] != tag)
					{
						vis[i1][j1] = tag;
						addFace(i1, j1, new Face(*newFace));
					}
				}
			}
			else
			{
				for (ch = getc(fp); ch !='\n' && ch != '\r'; ch = getc(fp));
				continue;
			}
		}

//		printf("X: %f %f  Z: %f %f\n", cityMinX, cityMaxX, cityMinZ, cityMaxZ);
//		cityMinX--;
//		cityMinZ--;	
		printf("initialize map border over\n");
	}

	void applyForce(const Vector2f &f, const Vector2f &r)
	{
		debug
		{
			printf("applyForce f :"); f.Print(0);
			printf(" r :"); r.Print();
		}

		combinedForce += f;
		static float someFactor = (sqr(carWidth)+sqr(carLength))*carMass/12;
		combinedAngularAcceleration -= det(f,r)/someFactor;
	}

	void calculateWheelForce(const Vector2f &wheelVector, const Vector2f &wheelDirection)
	{
		Vector2f tangentVelocity = wheelVector.Rotate(M_PI/2).Normalize() * carAngularVelocity;
		Vector2f relativeVelocity = wheelDirection*wheelVelocity - carVelocity - tangentVelocity;


		debug
		{
			printf("tangentVelocity: "); tangentVelocity.Print();
			printf("relativeVelocity: "); relativeVelocity.Print();
		}


		if (relativeVelocity.norm() < velocityEps)
		{
			if (abs(wheelVelocity) > 1e-3) applyForce(wheelDirection*carMass*frictionFactor, wheelVector);
		}
		else
		{
			Vector2f tmp = relativeVelocity.Normalize()*carMass*frictionFactor;
			applyForce(relativeVelocity.Normalize()*carMass*frictionFactor, wheelVector);
		}
		applyForce(-carVelocity*bounceFactor, wheelVector);
	}

	template <class T>
	void add(T &a, T b)
	{
		if (a*(a+b) < 0) a = 0;
		else a += b;
	}


    void Run()
    {
        GLUTBackendRun(this);
    }


	int sb = 0;

    virtual void RenderSceneCB()
	{

		debug
		{
			printf("\nrender\n");
		}

        m_scale += 0.01f;

        m_pGameCamera->OnRender();

		Vector3f toDirec3D = Vector3f(1.0f, 0.0f, 1.0f);
		toDirec3D.Rotate(car.yAngle, Vector3f(0.0f, 1.0f, 0.0f));
		toDirec3D.Normalize();
		Vector2f headDirection(toDirec3D.x, toDirec3D.z);
		Vector3f Right = Vector3f(0.0f, 1.0f, 0.0f).Cross(toDirec3D);
        Right.Normalize();
		Vector2f rightDirection(Right.x, Right.z);



		//wheelVelocity*************
		if (keyStatus['w'] == 1)
		{
			wheelVelocity += wheelAcceleration;
			if (wheelVelocity > maxWheelVelocity) wheelVelocity = maxWheelVelocity;
		}
		else
		{
			if (wheelVelocity > 0)
			{
				wheelVelocity -= wheelAcceleration;
				if (wheelVelocity < 0) wheelVelocity = 0;
			}
		}
		if (keyStatus['s'] == 1)
		{
			if (wheelVelocity > 0)
			{
				wheelVelocity -= wheelDeceleration;
				if (wheelVelocity < 0) wheelVelocity = 0;
			}
			else
			{
				wheelVelocity -= wheelDeceleration/2;
				if (wheelVelocity < -maxWheelReverseVelocity) wheelVelocity = -maxWheelReverseVelocity;
			}
		}
		else
		{
			if (wheelVelocity < 0)
			{
				wheelVelocity += wheelDeceleration;
				if (wheelVelocity > 0) wheelVelocity = 0;
			}
		}
		//**************************

		
		//wheelAngle****************
		if (keyStatus['a'] == 1)
		{
			wheelAngle -= turningVelocity;
			if (wheelAngle < -maxWheelAngle) wheelAngle = -maxWheelAngle;
		}
		else
		{
			if (wheelAngle < 0)
			{
				wheelAngle += turningVelocity * 2;
				if (wheelAngle > 0) wheelAngle = 0;
			}
		}
		if (keyStatus['d'] == 1)
		{
			wheelAngle += turningVelocity;
			if (wheelAngle > maxWheelAngle) wheelAngle = maxWheelAngle;
		}
		else
		{
			if (wheelAngle > 0)
			{
				wheelAngle -= turningVelocity * 2;
				if (wheelAngle < 0 ) wheelAngle = 0;
			}
		}
		//**************************

		Vector2f carCenter(car.center.x, car.center.z);
		Vector2f frontLeftWheelVector = -rightDirection*(carWidth*0.5) + headDirection*(carLength*0.42);
		Vector2f frontRightWheelVector = rightDirection*(carWidth*0.5) + headDirection*(carLength*0.42);
		Vector2f rearLeftWheelVector = -rightDirection*(carWidth*0.5) - headDirection*(carLength*0.35);
		Vector2f rearRightWheelVector = rightDirection*(carWidth*0.5) - headDirection*(carLength*0.35);
		Vector2f wheelDirection = headDirection.Rotate(wheelAngle);
		
		combinedForce = Vector2f(0,0);
		combinedAngularAcceleration = 0;
		calculateWheelForce(frontLeftWheelVector, wheelDirection);
		calculateWheelForce(frontRightWheelVector, wheelDirection);
		calculateWheelForce(rearLeftWheelVector, headDirection);
		calculateWheelForce(rearRightWheelVector, headDirection);
		
		debug
		{
			printf("wheelVelocity=%lf\n", wheelVelocity);
			printf("wheelAngle=%lf\n", wheelAngle);
			printf("combinedForce: "); combinedForce.Print();
			printf("conminefAngularAcceleration: %lf\n", combinedAngularAcceleration);
		}

		add(carVelocity.x, (combinedForce/carMass).x);
		add(carVelocity.y, (combinedForce/carMass).y);
		add(carAngularVelocity, combinedAngularAcceleration);
		add(wheelRoll, -wheelVelocity*27);
		if (carVelocity.norm() < 1e-3) carAngularVelocity = 0;

		debug
		{
			printf("carVelocity: "); carVelocity.Print();
			printf("carAngularVelocity: %lf\n", carAngularVelocity); 
		}

		Vector3f newCenter = Vector3f(car.center.x+carVelocity.x, car.center.y, car.center.z+carVelocity.y);

		float newYAngle = car.yAngle + ToDegree(carAngularVelocity);

		toDirec3D = Vector3f(1.0f, 0.0f, 1.0f);
		toDirec3D.Rotate(newYAngle, Vector3f(0.0f, 1.0f, 0.0f));
		toDirec3D.Normalize();
		Right = Vector3f(0.0f, 1.0f, 0.0f).Cross(toDirec3D);
        Right.Normalize();

		Vector3f fl, fr, rl, rr;
		fl = newCenter + toDirec3D*(carLength/2) - Right*(carWidth/2);
		fr = newCenter + toDirec3D*(carLength/2) + Right*(carWidth/2);
		rl = newCenter - toDirec3D*(carLength/2*0.9) - Right*(carWidth/2);
		rr = newCenter - toDirec3D*(carLength/2*0.9) + Right*(carWidth/2);

		if (inersectWithBarrier(fl.x,fl.z,fr.x,fr.z) || inersectWithBarrier(fl.x,fl.z,rl.x,rl.z) || inersectWithBarrier(fr.x,fr.z,rr.x,rr.z) || inersectWithBarrier(rl.x,rl.z,rr.x,rr.z))
		{
			carVelocity = Vector2f(0,0);
			carAngularVelocity = 0;
			wheelVelocity = 0;
		}
		else
		{
			car.yAngle = newYAngle;
			car.center = newCenter;
		}
//		printf("car.center: "); car.center.Print(); printf("\n");

		toDirec3D = Vector3f(1.0f, 0.0f, 1.0f);
		toDirec3D.Rotate(car.yAngle, Vector3f(0.0f, 1.0f, 0.0f));
		toDirec3D.Normalize();
		Right = Vector3f(0.0f, 1.0f, 0.0f).Cross(toDirec3D);
        Right.Normalize();

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		m_pEffect->Enable();

		const Vector3f Vaxis(0.0f, 1.0f, 0.0f);

		Vector3f camTarget = Vector3f(3.0f, -1.0f, 3.0f);
//		if (car.yAngle > 360.0f){ car.yAngle -= 360.0f; }

		viewYAngle += (car.yAngle-viewYAngle)*0.1;

		camTarget.Rotate(viewYAngle, Vaxis);

		
		camTarget.Rotate(cameraRotateH, Vaxis);

		Vector3f Haxis = Vaxis.Cross(camTarget);
		Haxis.Normalize();
		camTarget.Rotate(cameraRotateV, Haxis);
       

		camTarget.Normalize();

		Vector3f camPos = car.center - (camTarget * cameraDis);

		double splitAngle = angle(carVelocity.x, carVelocity.y, headDirection.x, headDirection.y);

//		camPos -= Right*splitAngle*carVelocity.norm()*5;

		checkIntersectionWithFace(car.center, camPos);

		SpotLight sl;
		sl.DiffuseIntensity = 0.9f;
		sl.Color = Vector3f(0.0f, 1.0f, 1.0f);
		sl.Position = car.center;
		sl.Direction = Vector3f(toDirec3D.x, 0.0, toDirec3D.z);
		sl.Attenuation.Linear = 0.1f;
		sl.Cutoff = 10.0f;

		m_pEffect->SetSpotLights(1, &sl);



		/*************************** Start ShadowMapPass ****************************/
		m_shadowMapFBO.BindForWriting();
		glClear(GL_DEPTH_BUFFER_BIT);

		m_ShadowMapEffect.Enable();

		Pipeline ps;
		ps.SetCamera(car.center + Vector3f(0.5, 0.0, 0.5), m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
		ps.SetOrthographicProj(m_shadowOrthoProjInfo);

		/* Draw City */
		ps.Scale(1.0f, 1.0f, 1.0f);
		ps.Rotate(0.0f, 0.0f, 0.0f);
		ps.WorldPos(0.0f, 0.0f, 0.0f);
		m_ShadowMapEffect.SetWVP(ps.GetWVOrthoPTrans());
		m_pMesh->Render();


		/* Draw Car */
		ps.Scale(carScale_s, carScale_s, carScale_s);
		ps.Rotate(0.0f, 45.0f - car.yAngle, 0.0f);
		ps.WorldPos(car.center-Vector3f(0.0f,0.35f,0.0f));
		m_ShadowMapEffect.SetWVP(ps.GetWVOrthoPTrans());
		if (displayCar) car_pMesh->Render();

		/* Draw rearWheel */
		ps.Scale(carScale_s, carScale_s, carScale_s);
		ps.Rotate(0.0f, 45.0f - car.yAngle, wheelRoll);
		ps.WorldPos(car.center - toDirec3D*carScale_s * 3700 + Vector3f(0.0f, carScale_s * 580, 0.0f));
		m_ShadowMapEffect.SetWVP(ps.GetWVOrthoPTrans());
		if (displayCar) rearWheel_pMesh->Render();

		/* Draw frontLeftWheel */
		ps.Scale(carScale_s, carScale_s, carScale_s);
		ps.Rotate(0.0f, 45.0f - car.yAngle - ToDegree(wheelAngle), wheelRoll);
		ps.WorldPos(car.center + toDirec3D*carScale_s * 5300 + Vector3f(0.0f, carScale_s * 580, 0.0f) - Right*carScale_s * 2700);
		m_ShadowMapEffect.SetWVP(ps.GetWVOrthoPTrans());
		if (displayCar) frontLeftWheel_pMesh->Render();

		/* Draw frontRightWheel */
		ps.Scale(carScale_s, carScale_s, carScale_s);
		ps.Rotate(0.0f, 45.0f - car.yAngle - ToDegree(wheelAngle), wheelRoll);
		ps.WorldPos(car.center + toDirec3D*carScale_s * 5300 + Vector3f(0.0f, carScale_s * 580, 0.0f) + Right*carScale_s * 2700);
		m_ShadowMapEffect.SetWVP(ps.GetWVOrthoPTrans());
		if (displayCar) frontRightWheel_pMesh->Render();



		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		/*************************** End**************************** */

		/*************************** Start RenderPass ***********************/
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		m_pEffect->Enable();
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_shadowMapFBO.BindForReading(SHADOW_TEXTURE_UNIT);


        Pipeline p;
        p.Scale(1.0f, 1.0f, 1.0f);
        p.Rotate(0.0f, 0.0f, 0.0f);
        p.WorldPos(0.0f, 0.0f, 0.0f);

						/********* Set LightSpace Position *************/
						p.SetOrthographicProj(m_shadowOrthoProjInfo);
						p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
						m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
						/***********************************************************/

		m_pGameCamera->SetCameraState(camPos, camTarget, Vector3f(0.0f, 1.0f, 0.0f));

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
        p.SetPerspectiveProj(m_persProjInfo);
        m_pEffect->SetWVP(p.GetWVPTrans());
        m_pEffect->SetWorldMatrix(p.GetWorldTrans());
        m_pEffect->SetDirectionalLight(m_directionalLight);
        m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
        m_pEffect->SetMatSpecularIntensity(0.0f);
        m_pEffect->SetMatSpecularPower(0);

        m_pMesh->Render();




		/* car *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle, 0.0f);
		p.WorldPos(car.center-Vector3f(0.0f,0.35f*carScale/0.0007f,0.0f));

					/********* Set LightSpace Position *************/
					p.SetOrthographicProj(m_shadowOrthoProjInfo);
					p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
					m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
					/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) car_pMesh->Render();
		/*********************/


		/* rearWheel *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle, wheelRoll);
		p.WorldPos(car.center - toDirec3D*carScale*3700 + Vector3f(0.0f, carScale*580, 0.0f));

					/********* Set LightSpace Position *************/
					p.SetOrthographicProj(m_shadowOrthoProjInfo);
					p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
					m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
					/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) rearWheel_pMesh->Render();
		/*********************/
		

		/* frontLeftWheel *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle-ToDegree(wheelAngle*2), wheelRoll);
		p.WorldPos(car.center + toDirec3D*carScale*5300 + Vector3f(0.0f, carScale*580, 0.0f) - Right*carScale*2700);

					/********* Set LightSpace Position *************/
					p.SetOrthographicProj(m_shadowOrthoProjInfo);
					p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
					m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
					/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) frontLeftWheel_pMesh->Render();
		/*********************/

		/* frontRightWheel *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle-ToDegree(wheelAngle*2), wheelRoll);
		p.WorldPos(car.center + toDirec3D*carScale*5300 + Vector3f(0.0f, carScale*580, 0.0f) + Right*carScale*2700);

						/********* Set LightSpace Position *************/
						p.SetOrthographicProj(m_shadowOrthoProjInfo);
						p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
						m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
						/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) frontRightWheel_pMesh->Render();
		/*********************/


		/* frontLeftBrake *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle-ToDegree(wheelAngle*2), 0.0f);
		p.WorldPos(car.center + toDirec3D*carScale*5300 + Vector3f(0.0f, carScale*580, 0.0f) - Right*carScale*2700);

						/********* Set LightSpace Position *************/
						p.SetOrthographicProj(m_shadowOrthoProjInfo);
						p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
						m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
						/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) frontLeftBrake_pMesh->Render();
		/*********************/

		/* frontRightBrake *******************/
		p.Scale(carScale, carScale, carScale);
		p.Rotate(0.0f, 45.0f-car.yAngle-ToDegree(wheelAngle*2), 0.0f);
		p.WorldPos(car.center + toDirec3D*carScale*5300 + Vector3f(0.0f, carScale*580, 0.0f) + Right*carScale*2700);

							/********* Set LightSpace Position *************/
							p.SetOrthographicProj(m_shadowOrthoProjInfo);
							p.SetCamera(car.center, m_directionalLight.Direction, Vector3f(0.0f, 1.0f, 0.0f));
							m_pEffect->SetLightWVP(p.GetWVOrthoPTrans());
							/***********************************************************/

		p.SetCamera(m_pGameCamera->GetPos(), m_pGameCamera->GetTarget(), m_pGameCamera->GetUp());
		p.SetPerspectiveProj(m_persProjInfo);
		m_pEffect->SetWVP(p.GetWVPTrans());
		m_pEffect->SetWorldMatrix(p.GetWorldTrans());
		m_pEffect->SetDirectionalLight(m_directionalLight);
		m_pEffect->SetEyeWorldPos(m_pGameCamera->GetPos());
		m_pEffect->SetMatSpecularIntensity(0.0f);
		m_pEffect->SetMatSpecularPower(0);
		if (displayCar) frontRightBrake_pMesh->Render();
		/*********************/
			

		if (displaySkyBox) m_pSkyBox->Render();

        glutSwapBuffers();
    }

    void KeyboardCB(OGLDEV_KEY OgldevKey, OGLDEV_KEY_STATE OgldevKeyState = OGLDEV_KEY_STATE_PRESS)
    {
        switch (OgldevKey) {
		case OGLDEV_KEY_q:
		if (carScale > 0.0006f)
		{
			carScale = 0.0001f;
			carScale_s = 0.0001f;
			carLength /= 7;
			carWidth /= 7;
		}
		else
		{
			carScale = 0.0007f;
			carScale_s = 0.0008f;
			carLength *= 7;
			carWidth *= 7;
			maxWheelAngle -= 0.2;
			frictionFactor /= 3;
		}
		break;
        default:
            m_pGameCamera->OnKeyboard(OgldevKey);				
        }
    }

    virtual void PassiveMouseCB(int x, int y)
    {
		cameraRotateH = (x-WINDOW_WIDTH/2.0)/WINDOW_WIDTH*360;
		if (y < WINDOW_HEIGHT/2)
		{
			y = (0.37+(y*1.0/WINDOW_HEIGHT)*0.26)*WINDOW_HEIGHT;
		}
		cameraRotateV = (y-WINDOW_HEIGHT/2.0)/WINDOW_HEIGHT*80;
    }

	void MouseCB(OGLDEV_MOUSE OgldevMouse, OGLDEV_KEY_STATE OgldevKeyState, int x, int y)
	{
	}

	void MouseWheelCB(int direction)
	{
		if (direction == 1)
		{
			if (cameraDis-1 > 6) cameraDis -= 1;
		}
		else
		{
			if (cameraDis+1 < 25) cameraDis += 1;
		}
	}

private:

    LightingTechnique* m_pEffect;
    Camera* m_pGameCamera;
    float m_scale;
    DirectionalLight m_directionalLight;
    Mesh* m_pMesh, *car_pMesh, *rearWheel_pMesh, *frontLeftWheel_pMesh, *frontRightWheel_pMesh, *frontLeftBrake_pMesh, *frontRightBrake_pMesh;
    PersProjInfo m_persProjInfo;	
	tankModel car;
	SkyBox* m_pSkyBox;

	ShadowMapFBO m_shadowMapFBO;
	OrthoProjInfo m_shadowOrthoProjInfo;
	ShadowMapTechnique m_ShadowMapEffect;
	SpotLight m_spotlight;
};


int main(int argc, char** argv)
{
	printf("start main\n");
    Magick::InitializeMagick(*argv);
    GLUTBackendInit(argc, argv, true, false);


    if (!GLUTBackendCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, false, "Tutorial 22")) {
        return 1;
    }

    Tutorial22* pApp = new Tutorial22();

    if (!pApp->Init()) {
        return 1;
    }

    pApp->Run();

    delete pApp;
 
    return 0;
}
