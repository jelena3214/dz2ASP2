#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <fstream>

using namespace std;

const int nn = 3; //red stabla br pokazivaca
const int order = 2; //2 * floor((2 * nn - 2) / 3) + 1; OVO JE MAX BR STRINGOVA, SAMO KOREN IMA TOLIKO OSTALI
//SE POPUNJAVAJU DO NN-1

struct Node {

private:

	void init() {
		for (int i = 0; i < order; i++) {
			data[i] = nullptr;
			pointers[i] = nullptr;
		}
		pointers[order] = nullptr;
		currElems = 0;
		father = nullptr;
	}

public:
	string* data[order];
	Node* pointers[order + 1];
	int leaf;
	int root;
	int currElems;
	Node* father;
	Node() {
		init();
	}

	/*Node(const Node& n);
	Node(Node&& n);
	Node& operator=(const Node& n); {
		leaf = n.leaf;
		root = n.root;
		currElems = n.currElems;
		father = n.father;
		for (int i = 0; i < order; i++) {
			data[i] = new string{ *n.data[i] };
			pointers[i] = n.pointers[i];
		}
		pointers[order] = n.pointers[order];
	}
	//Node& operator=(Node&& n);*/
	~Node() {
		for (int i = 0; i < order; i++) {
			delete data[i];
			delete pointers[i];
		}
		delete pointers[order];
	}
	int maxPointers() {
		if (root)
			return 2 * floor((2 * nn - 2) / 3) + 1;
		else return nn;
	}
	int minPointers() {
		if (root) return 2;
		else return ceil((2 * nn - 1) / 3);
	}
	int maxKeys() { return maxPointers() - 1; }

	bool isFull() { return currElems == maxPointers() - 1; }

	friend ostream& operator<<(ostream& os, Node* root) {
		queue<Node*>qu;
		qu.push(root);

		Node* gran = new Node;
		gran->data[0] = new string{ "XXXXX" };
		qu.push(gran);

		while (!qu.empty()) {
			Node* tmp = qu.front();
			qu.pop();

			if (qu.size() == 0) {
				return os;
			}
			if (*tmp->data[0] == "XXXXX") {
				os << endl;
				qu.push(gran);
			}
			else {
				os << "[";
				for (int i = 0; i < tmp->currElems; i++) {
					os << *tmp->data[i] << " ";
				}
				os << "]   ";
				
				for (int i = 0; i <= tmp->currElems; i++) {
					if (tmp->pointers[i])qu.push(tmp->pointers[i]);
				}
				
			}
		}
		return os;
	}
};
void insertNode(Node* root, string d);

Node* rightBrother(Node* tmp, int& pointerIndex, int& keyIndex) {
	Node* fath = tmp->father;
	if (fath == nullptr)return fath;
	for (int i = 0; i <= fath->currElems; i++) {
		if (fath->pointers[i] == tmp)pointerIndex = i;
	}
	keyIndex = pointerIndex - 1;
	if (pointerIndex + 1 <= fath->currElems)return fath->pointers[pointerIndex + 1];
	return nullptr;
}

Node* leftBrother(Node* tmp, int& pointerIndex, int& keyIndex) {
	Node* fath = tmp->father;
	if (fath == nullptr)return fath;
	for (int i = 0; i <= fath->currElems; i++) {
		if (fath->pointers[i] == tmp)pointerIndex = i;
	}
	keyIndex = pointerIndex - 1;
	if (pointerIndex - 1 >= 0)return fath->pointers[pointerIndex - 1];
	return nullptr;
}

int smallerKeys(Node* root, string key) {
	queue<Node*>qu;
	qu.push(root);
	vector<string>keys;
	int num = 0, i = 0;
	while (!qu.empty()) {
		Node* tmp = qu.front();
		qu.pop();
		
		for (int i = 0; i < tmp->currElems; i++) {
			keys.push_back(*tmp->data[i]);
		}

		for (int i = 0; i <= tmp->currElems; i++) {
			if (tmp->pointers[i])qu.push(tmp->pointers[i]);
		}
	}
	sort(keys.begin(), keys.end());
	while (keys[i++] < key)num++;
	return num;
}

Node* searchKey(Node* root, string key, int& pos) {
	Node* temp = root, * find = root;
	int greater = 0, flag = 1;//flag znaci da je veci od svih
	while (temp) {
		find = temp;
		for (int i = 0; i < temp->currElems; i++) {
			if (*temp->data[i] > key) {
				flag = 0;
				temp = temp->pointers[i];
				if (!temp)return nullptr;
				break;
			}
			else if(*temp->data[i] == key){
				pos = i;
				return temp;
			}
			else {
				greater++;
			}
		}
		if (greater == temp->currElems && flag)temp = temp->pointers[temp->currElems]; //ako je veci od svih
		greater = 0;
		flag = 1;
	}
	return nullptr;
}

Node* findLeafToInsert(Node* root, string d, int& pos) {
	Node* temp = root, * find = root;
	pos = -1;
	int greater = 0, flag = 1;//flag znaci da je veci od svih
	while (temp) {
		find = temp;
		for (int i = 0; i < temp->currElems; i++) {
			if (*temp->data[i] > d) {
				flag = 0;
				temp = temp->pointers[i], pos = i;
				if (!temp)return find;
				break;
			}
			else {
				greater++;
			}
		}
		if (greater == temp->currElems && flag)pos = temp->currElems, temp = temp->pointers[temp->currElems]; //ako je veci od svih
		greater = 0;
		flag = 1;
	}
	return find;
}

void deleteKeys(Node* nod) {
	for (int i = 0; i < nod->currElems; i++) {
		delete nod->data[i];
		nod->data[i] = nullptr;
	}
	//nod->currElems = 0;
}

int numOfPointers(Node* nod) {
	Node* tmp = nod;
	int i = 0;
	while (tmp->pointers[i++]);
	return i;
}

bool greaterArray(string** s1, string** s2, int len1, int len2) {
	/*if (len1 == 1 && len2 == 1) {
		if (*s1[0] > *s2[0])return true;
		else return false;
	}
	else {
		int first = 0, second = 0;
		for (int i = 0; i < len1; i++) {
			if (*s1[i] > *s2[i])first++;
			else second++;
		}
		if (first == second)return true;
		else if (first > second)return true;
		else return false;
	}*/
	if (*s1[0] > *s2[0])return true;
	else return false;
}

void sortVect(vector<Node*>&point) {
	for (int i = 0; i < point.size()-1; i++) {
		for (int j = i + 1; j < point.size(); j++) {
			if (greaterArray(point[i]->data, point[j]->data, point[i]->currElems, point[j]->currElems)) {//da li je prvi niz veci od drugog
				Node* tmp = point[i];
				point[i] = point[j];
				point[j] = tmp;
			}
		}
	}
}

vector<vector<Node*>> levelOrder(Node* root) {
		vector < vector <Node*> > ans;
		if (!root)return ans;
		queue <Node*> q;
		q.push(root);
		while (!q.empty()) {
			int sz = q.size();
			vector<Node*> temp;
			while (sz--) {
				Node* curr = q.front();
				temp.push_back(curr);
				q.pop();
				for (int i = 0; i <= curr->currElems; i++) {
					if(curr->pointers[i])q.push(curr->pointers[i]);
				}
			}
			ans.push_back(temp);
		}
		return ans;
}

void nodeSeparating(Node* separate, string d, int pos, int rl, Node* root) {//pos nam pozicija prelomnog cvora, rl = 1 kad ima desnog brata
	int n = separate->currElems + (separate->root ? 1 : 2);
	int inserted = 0;

	int firstNum = floor((2 * nn - 2) / 3);
	int secondNum = floor((2 * nn - 1) / 3);
	int thirdNum = floor((2 * nn) / 3);

	vector<string> allData;

	if (separate->root) {//koren lomimo na 2 cvora
		for (int i = 0; i < separate->currElems; i++)allData.push_back(*separate->data[i]);
		allData.push_back(d);
		sort(allData.begin(), allData.end());

		deleteKeys(separate);
		int mid = allData.size() / 2;
		Node* firstNode = new Node;
		for (int i = 0; i < mid; i++) {
			firstNode->data[i] = new string{ allData[i] };
			firstNode->currElems++;
		}
		separate->data[0] = new string{ allData[mid] };
		separate->pointers[0] = firstNode;
		separate->pointers[0]->leaf = 1;
		separate->pointers[0]->root = 0;
		firstNode->father = separate;
		separate->currElems -= firstNode->currElems - 1;

		Node* secondNode = new Node;
		for (int i = 0; i < allData.size() - mid - 1; i++) {
			secondNode->data[i] = new string{ allData[i + mid + 1] };
			secondNode->currElems++;
		}
		secondNode->father = separate;
		separate->pointers[1] = secondNode;
		separate->pointers[1]->leaf = 1;
		separate->pointers[1]->root = 0;
		separate->currElems -= secondNode->currElems;
	}
	else {//PRELAMANJE CVORA KOJI NIJE KOREN?
		int brotherPos = (rl == 1 ? pos + 1 : pos - 1); //pos je pzicija pointera 
		Node* fath = separate->father;
		Node* brother = fath->pointers[brotherPos];
		for (int i = 0; i < separate->currElems; i++)allData.push_back(*separate->data[i]);
		for (int i = 0; i < brother->currElems; i++)allData.push_back(*brother->data[i]);
		allData.push_back(d);
		int keyInd;// = (pos == 0) ? 0 : pos - 1;
		if (rl == 1) {
			keyInd = pos;
		}
		else {
			keyInd = pos - 1;
			if (keyInd < 0)keyInd = 0;
		}
		allData.push_back(*fath->data[keyInd]);
		sort(allData.begin(), allData.end());

		deleteKeys(separate);
		deleteKeys(brother);
		separate->currElems = 0;
		brother->currElems = 0;
		if (rl) {
			for (int i = 0; i < firstNum; i++)separate->data[i] = new string{ allData[i] }, separate->currElems++;
			for (int i = firstNum + 1, k = 0; i < secondNum + firstNum + 1; i++)brother->data[k++] = new string{ allData[i] }, brother->currElems++;
		}
		else {
			for (int i = 0; i < firstNum; i++)brother->data[i] = new string{ allData[i] }, brother->currElems++;
			for (int i = firstNum + 1, k = 0; i < secondNum + firstNum + 1; i++)separate->data[k++] = new string{ allData[i] }, separate->currElems++;
		}

		Node* thirdNode = new Node;
		for (int i = secondNum + 2 + firstNum, k = 0; i < allData.size(); i++)thirdNode->data[k++] = new string{ allData[i] }, thirdNode->currElems++;
		thirdNode->father = fath;
		thirdNode->root = 0;
		thirdNode->leaf = 1;

		vector<string> nodeData;
		nodeData.push_back(allData[firstNum]);
		nodeData.push_back(allData[firstNum + secondNum + 1]);
		sort(nodeData.begin(), nodeData.end());
		delete fath->data[keyInd];
		fath->data[keyInd] = nullptr;
		if (fath->currElems + 1 <= fath->maxKeys()) {
			fath->data[keyInd] = new string{ nodeData[0] };
			//sve pokazivace posle pos za 1 gore i tako i kljuceve posle keyind
			for (int i = fath->currElems; i > pos; i--) fath->pointers[i + 1] = fath->pointers[i];
			for (int i = fath->currElems-1; i > keyInd+1; i--) fath->data[i + 1] = fath->data[i];
			fath->data[keyInd+1] = new string{ nodeData[1] };
			fath->pointers[keyInd + 2] = thirdNode;
			fath->currElems++;
		}
		else {
			vector<Node*>rootpoint;
			Node* rootPointers[order + 2]; //povecaj ovo
			for (int i = 0; i <= fath->currElems; i++)rootPointers[i] = fath->pointers[i],rootpoint.push_back(fath->pointers[i]), fath->pointers[i] = nullptr;//preruzimamo njegove sinove
			rootPointers[fath->currElems + 1] = thirdNode;
			rootpoint.push_back(thirdNode);
			//moramo da preuzmemo sve od njegove brace od leva pa na desno:

			vector<vector<Node*>>vec = levelOrder(root);
			Node* tmp = fath;
			int lev = 0;
			while (tmp->root != 1) {
				tmp = tmp->father;
				lev++;
			}
			if (vec.size() == 1)lev = 0;
			else lev++;
			vector<Node*>v = vec[lev];
			if (lev == 0)v.clear();
			for (int i = 0; i < rootpoint.size(); i++)v.push_back(rootpoint[i]);
			sortVect(v);

			rootpoint = v;

			fath->currElems--;
			//pravimo 2 nova cvora, prelamamo koren
			insertNode(fath, nodeData[0]);
			insertNode(fath, nodeData[1]);
			Node* left, * right;
			if (fath->root == 1) {
				int k = 0;
				for (int i = 0; i <= fath->currElems; i++) {
					for (int j = 0; j <= fath->pointers[i]->currElems; j++) {
						fath->pointers[i]->pointers[j] = rootpoint[k++];
						fath->pointers[i]->pointers[j]->father = fath->pointers[i];
					}
				}
			}
			else {
				Node* tmp = fath->father;
				int i = 0; 
				vector<Node*>brothers;

				int lev = 1;
				while (tmp->root !=1) {
					tmp = tmp->father;
					lev++;
				}
				//treba ubaciti sve iz nivoa u brothers
				vector<vector<Node*>>vec = levelOrder(root);
				brothers = vec[lev];



				int k = 0, l = 0;
				while(brothers.size()>l) {
					for (int j = 0; j <= brothers[l]->currElems; j++) {
						brothers[l]->pointers[j] = rootpoint[k++];
						brothers[l]->pointers[j]->father = brothers[l];
					}
					l++;
				}
			}
		}

	}

}

void nodePouring(Node* pour, Node* brother, int pos, string d, int rl) {//pos lokacija prelomnog
	vector<string> allData;
	Node* fath = pour->father;
	for (int i = 0; i < pour->currElems; i++)allData.push_back(*pour->data[i]);
	for (int i = 0; i < brother->currElems; i++)allData.push_back(*brother->data[i]);
	int pozition; // = (rl == 1 ? pos : pos - 1);
	if (rl == 1) {
		pozition = pos;
	}
	else {
		pozition = pos - 1;
		if (pozition < 0)pozition = 0;
	}
	
	string s = *pour->father->data[pozition];
	allData.push_back(*pour->father->data[pozition]);
	allData.push_back(d);

	sort(allData.begin(), allData.end());

	int mid = (allData.size()-1) / 2; //zbog parnih brojeva da zaokruzi na donji
	delete fath->data[pozition];
	fath->data[pozition] = new string{ allData[mid] };

	deleteKeys(pour);
	deleteKeys(brother);
	pour->currElems = 0;
	brother->currElems = 0;

	if (rl) {
		for (int i = 0; i < mid; i++)pour->data[i] = new string{ allData[i] }, pour->currElems++;
		for (int i = mid + 1, k = 0; i < allData.size(); i++)brother->data[k++] = new string{ allData[i] }, brother->currElems++;
	}
	else {
		for (int i = 0; i < mid; i++)brother->data[i] = new string{ allData[i] }, brother->currElems++;
		for (int i = mid + 1, k = 0; i < allData.size(); i++)pour->data[k++] = new string{ allData[i] }, pour->currElems++;
	}

}

void insertNode(Node* root, string d) {

	int positionOfSon = 0;

	if (root->data[0] == nullptr) {//prvi string umecemo
		root->data[0] = new string{ d };
		root->currElems += 1;
	}
	else {
		Node* place = findLeafToInsert(root, d, positionOfSon);
		if (place->currElems < place->maxKeys()) {//ima mesta u cvoru
			for (int j = place->currElems; j > positionOfSon; j--)place->data[j] = place->data[j - 1];
			place->data[positionOfSon] = new string{ d };
			place->currElems++;
		}
		else {//prelivanje ako nema bracu prelamanje
			int flag = 0, rl = 0, separatingrl = 0, pointerIndex = 0, keyIndex = 0;
			if (place->root) {//kod njega nema prelivanja samo prelamanje
				nodeSeparating(place, d, positionOfSon,0, root);
			}
			else {
				Node* help = nullptr;
				//prelivanje ako moze, ako ne prelamamo
				int hasright = 0;
				if (rightBrother(place, pointerIndex, keyIndex)) {//da li je u opsegu
					separatingrl = 1;
					if (!rightBrother(place, pointerIndex, keyIndex)->isFull()) {
						help = rightBrother(place, pointerIndex, keyIndex);
						hasright = 1; //za prelamanje prioritet desni kao ima oba
						flag = 1;
						rl = 1; 
					}
				}
				if (leftBrother(place, pointerIndex, keyIndex)) { //da li je u opsegu
					//separatingrl = 0;
					if (!leftBrother(place, pointerIndex, keyIndex)->isFull()&&!hasright) {
						help = leftBrother(place, pointerIndex, keyIndex);
						flag = 1;
					}
				}
				if (flag) {
					nodePouring(place, help, pointerIndex, d, rl);
				}
				else {
					nodeSeparating(place, d, pointerIndex, separatingrl, root);
				}
			}

		}
	}
}

bool isLeaf(Node* tmp) {
	int l = 0;
	for (int i = 0; i <= tmp->currElems; i++) {
		if (tmp->pointers[i] != nullptr)l++;
	}
	return (l == 0);
}

void loan(Node* curr, Node* brother, int keyInd) {//pozicija razdvojnog kljuca
	Node* fath = curr->father;
	//brisemo mu taj sto treba i ubacujemo pozajmicu
	int i = 0;
	if (curr->data[0] == nullptr)i = 0;
	else{ while (fath->data[keyInd] > curr->data[i++]); }//ide na kraj
	curr->data[i] = fath->data[keyInd];
	fath->data[keyInd] = brother->data[0];
	for (int i = 0; i < brother->currElems; i++)brother->data[i] = brother->data[i + 1];
	brother->currElems--;
}

void loanLeft(Node* curr, Node* brother, int keyInd) {//pozicija razdvojnog kljuca
	Node* fath = curr->father;
	//brisemo mu taj sto treba i ubacujemo pozajmicu
	int i = 0;
	if (curr->data[0] == nullptr)i = 0;
	else { //stavljamo ga na pocetak najveci je
		for (int i = curr->currElems - 1; i >= 0; i--)curr->data[i + 1] = curr->data[i]; //shift right
	}
	curr->data[i] = fath->data[keyInd];
	fath->data[keyInd] = brother->data[brother->currElems-1];
	brother->data[brother->currElems - 1] = nullptr;
	brother->currElems--;
}

bool helpFromBrother(Node* curr) {
	int pInd = 0, kInd = 0;
	Node* rightB = rightBrother(curr, pInd, kInd);
	Node* leftB = leftBrother(curr, pInd, kInd);
	if (kInd < 0)kInd = 0;
	if (rightB) { //ako postoji
		if (rightB->currElems - 1 > rightB->minPointers() - 1) {//ima doiovljno
			if (pInd == rightB->currElems)kInd = pInd - 1;
			else kInd = pInd;
			loan(curr, rightB, kInd);
			return true;
		}

	}
	else if (leftB) {
		if (leftB->currElems - 1 > leftB->minPointers() - 1) {
			if (pInd == leftB->currElems)kInd = pInd - 1;
			loanLeft(curr, leftB, kInd);
			return true;
		}
	}
	Node* right2B = nullptr, * left2B = nullptr;
	if(rightB)right2B = rightBrother(rightB, pInd, kInd);
	if(left2B)left2B = leftBrother(rightB, pInd, kInd);//ZA NULL
	if (rightB && right2B) { //ako postoji
		if (right2B->currElems - 1 > right2B->minPointers() - 1) {//ima doiovljno
			if (pInd == right2B->currElems)kInd = pInd - 1;
			else kInd = pInd;

			Node* fath = right2B->father;
			string* tmp = fath->data[kInd];
			fath->data[kInd] = right2B->data[0];
			right2B->data[0] = nullptr;
			for (int i = 0; i < right2B->currElems; i++)right2B->data[i] = right2B->data[i + 1];
			right2B->currElems--;
			rightB->data[rightB->currElems] = tmp;
			rightB->currElems++;
			tmp = fath->data[kInd - 1];
			fath->data[kInd - 1] = rightB->data[0];
			for (int i = 0; i < rightB->currElems; i++)rightB->data[i] = rightB->data[i + 1];
			rightB->currElems--;
			curr->data[curr->currElems-1] = tmp;
			//curr->currElems++;
			return true;
		}

	}
	else if (leftB && left2B) {
		if (left2B->currElems - 1 > left2B->minPointers() - 1) {
			if (pInd == left2B->currElems)kInd = pInd - 1;
			Node* fath = left2B->father;
			string* tmp = fath->data[kInd];
			fath->data[kInd] = left2B->data[left2B->currElems - 1];
			left2B->data[left2B->currElems - 1] = nullptr;
			left2B->currElems--;
				
			for (int i = leftB ->currElems - 1; i >= 0; i--)leftB->data[i + 1] = leftB->data[i]; //SHIFT RIGHT
			leftB->data[0] = tmp;
			leftB->currElems++;

			tmp = fath->data[kInd+1];
			fath->data[kInd+1] = leftB->data[leftB->currElems - 1];
			leftB->data[leftB->currElems - 1] = nullptr;
			leftB->currElems--;

			for (int i = curr->currElems - 1; i >= 0; i--)curr->data[i + 1] = curr->data[i]; //SHIFT RIGHT
			curr->data[0] = tmp;
			//curr->currElems++;

			return true;
		}
	}return false;
}

void merging(Node* curr, int p, Node*& root) {

	//trazimo 2 brata;
	int leftflag = 0, rightflag = 0, rlflag = 0;
	Node* brother1  = nullptr, * brother2 = nullptr;
	int ind1, ind2;
	int pIndr = 0, kIndr = 0, pIndl = 0, kIndl = 0;
	Node* rightB = rightBrother(curr, pIndr, kIndr);
	Node* leftB = leftBrother(curr, pIndl, kIndl);
//PROVERI OVE PIND KIND ITD	
	if (rightB) {
		if (leftB) {
			brother1 = rightB;
			brother2 = leftB;
			ind1 = pIndr;
			ind2 = pIndl;
			rlflag = 1;
		}
		else {
			int pIndr1 = 0, kIndr1 = 0;
			Node* right2B = rightBrother(rightB, pIndr1, kIndr1);
			if (right2B) {
				brother1 = rightB;
				brother2 = right2B;
				ind1 = pIndr;
				ind2 = pIndr1;
				rightflag = 1;
			}
		}
	}
	else if (leftB) {
		int pIndl1 = 0, kIndl1 = 0;
		Node* left2B = leftBrother(leftB, pIndl1, kIndl1);
		if (left2B) {
			brother1 = leftB;
			brother2 = left2B;
			ind1 = kIndl;
			ind2 = kIndl1;
			leftflag = 1;
		}
	}
	if(rlflag == 0 && rightflag == 0 && leftflag == 0) {
		curr->data[p] = nullptr;
		vector<Node*>allNodes;
		Node* tmp = curr;
		while (tmp->root != 1) {
			if (helpFromBrother(tmp))break;
			else {
				Node* fath = tmp->father;
				vector<Node*>vec;
				vector<string>allData;
				int pos = 0, pos1 = 0;
				Node* leftB = leftBrother(tmp, pos1, pos);
				Node* rightB = rightBrother(tmp, pos1, pos);
				Node* br = (leftB == nullptr ? rightB : leftB);

				for (int i = 0; i <= fath->currElems; i++) {
					if(fath->pointers[i])vec.push_back(fath->pointers[i]);
				}
				for (int i = 0; i < fath->currElems; i++) {
					if (fath->data[i])allData.push_back(*fath->data[i]);
					fath->data[i] = nullptr;
				}
				for (int i = 0; i < vec.size(); i++) {
					for (int j = 0; j < vec[i]->currElems; j++) {
						if(vec[i]->data[j])allData.push_back(*vec[i]->data[j]), fath->currElems--;
						
					}
				}
				sort(allData.begin(), allData.end());
				
				deleteKeys(fath);
				fath->currElems = 0;
				br->currElems = 0;
				deleteKeys(br);
				//stavljamo u brata kljuceve ostale
				for (int i = 0; i < allData.size(); i++) {
					br->data[i] = new string{ allData[i] };
					br->currElems++;
				}
				allNodes.push_back(br);
			}
			tmp = tmp->father;
			//ef bd poubaci a i c i spoji
		}
		allNodes[1]->pointers[2] = allNodes[0];
		for (int i = 0; i <= allNodes[0]->currElems; i++)allNodes[0]->pointers[i] = nullptr;
		allNodes[0]->father = allNodes[1];
		allNodes[1]->root = 1;
		allNodes[1]->leaf = 0;
		allNodes[1]->father = nullptr;
		root = allNodes[1];
		/*for (int i = 0; i < allNodes[1]->currElems; i++) {
			root->pointers[i] = allNodes[1]->pointers[i];
			root->data[i] = allNodes[1]->data[i];
			root->pointers[i]->currElems = allNodes[1]->pointers[i]->currElems;
		}
		root->pointers[allNodes[1]->currElems] = allNodes[1]->pointers[allNodes[1]->currElems];
		root->pointers[allNodes[1]->currElems]->currElems = allNodes[1]->pointers[allNodes[1]->currElems]->currElems;
		root->currElems = allNodes[1]->currElems;*/
		return;
	}

	vector<string>allData;
	Node* fath = curr->father;
	int index = 0; //pozicija curr-a u ocu
	for (int i = 0; i <= fath->currElems; i++)if (fath->pointers[i] == curr)index = i;

	for (int i = 0; i < brother1->currElems; i++)allData.push_back(*brother1->data[i]);
	for (int i = 0; i < brother2->currElems; i++)allData.push_back(*brother2->data[i]);
	allData.push_back(*fath->data[ind1]);
	allData.push_back(*fath->data[ind2]);
	fath->currElems -= 2;
	fath->data[ind1] = fath->data[ind2] = nullptr;

	//REMOVE NULLPTR IZ NIZA

	sort(allData.begin(), allData.end());
	int mid = (allData.size() - 1) / 2;
	int pos = (fath->currElems == 0 ? 0 : fath->currElems - 1);
	fath->data[pos] = new string{ allData[mid] };
	fath->currElems++;
	//jedan cvor brisemo
	deleteKeys(brother1);
	deleteKeys(brother2);
	brother1->currElems = brother2->currElems = 0;
	int k = 0;
	for (int i = 0; i < mid; i++)brother1->data[i] = new string{ allData[k++] }, brother1->currElems++;
	for (int i = k + 1, l = 0; i < allData.size(); i++) {
		brother2->data[l++] = new string{ allData[i] };
		brother2->currElems++;
	}
	brother1->father = fath;
	brother2->father = fath;
	int pozition = (ind1 < ind2 ? ind1 : ind2);
	int pozition2 = (ind1 < ind2 ? ind2 : ind1);
	fath->pointers[pozition] = brother1;
	fath->pointers[pozition2] = brother2;
	
	fath->pointers[index] = nullptr;
	//REMOVE NULLPTR IZ NIZA
	//AKO SE DESI DA U OCU SAD IMA MINIMALNO REKURZIJA;

}

Node* findSuccesor(Node* place, int pos) {
	Node* cur = place->pointers[pos + 1];
	while (!isLeaf(cur))
		cur = cur->pointers[0];

	return cur;
}

void deleteNode(Node*& root, string del) {
	//brisemo iz lista
	//brisemo iz cvora, premestamo u list;
	int pos = 0;
	Node* place = searchKey(root, del, pos); //ako je nullptr
	int index = 0;
	for (int i = 0; i < place->currElems; i++) {//trazimo poziciju
		if (*place->data[i] == del)index = i;
	}
	if (isLeaf(place)) {//onda vidimo da li moze da se samo ukloni
		if (place->currElems - 1 > place->minPointers() - 1) {//moze samo da se ukloni
			
			for (int i = index; i < place->currElems; i++)place->data[i] = place->data[i + 1]; //shift na levo
			place->currElems--;
		}
		else {
			//place->data[index] = nullptr;
			if(helpFromBrother(place))return;
			merging(place, index, root);

		}
	}
	else {//dovodimo ga u list
		Node* succ = findSuccesor(place, pos);
		//Node* rightP = place->pointers[pos+1];
		string* tmp = place->data[index];
		place->data[index] = succ->data[0];
		succ->data[0] = tmp;
		if (succ->currElems - 1 > succ->minPointers() - 1) {//moze samo da se ukloni
			for (int i = 0; i < succ->currElems; i++)succ->data[i] = succ->data[i + 1]; //shift na levo
			succ->currElems--;
		}
		else {
			if (helpFromBrother(succ))return;
			merging(succ, index, root);
		}
		
	}

}

//DATOTEKA
void data(Node* root) {
	ifstream MyReadFile;
	MyReadFile.open("t.txt");

	string key, line;
	vector<string> translations;
	struct TreeNode* node;

	while (getline(MyReadFile, key)) {
		getline(MyReadFile, line);
		insertNode(root, line);
	}

	MyReadFile.close();
}

int main() {
	Node* root = new Node;
	root->leaf = 0;
	root->root = 1;
	while (1) {
		cout << "===================================================== \n";
		cout << " \t\tMENI \t \n ";
		cout << "===================================================== \n";
		cout << " 1.Kreiraj stablo\n";
		cout << " 2.Ubaci element\n";
		cout << " 3.Pronadji kljuc \n";
		cout << " 4.Brisi kljuc \n";
		cout << " 5.Svi manji kljucevi od zadatog \n";
		cout << " 6.Ispisi stablo \n";
		cout << " 7.Kraj programa \n";
		int choice;
		cin >> choice;
		
		if (choice == 1) {
			root->leaf = 0;
			root->root = 1;
		}
		if (choice == 2) {
			string tmp;
			cout << " Unesi string: ";
			cin >> tmp;
			insertNode(root, tmp);
		}
		if (choice == 3) {
			string tmp;
			cout << " Unesi string: ";
			cin >> tmp;
			int pos;
			if (searchKey(root, tmp, pos)) {
				cout << "Postoji u stalu\n";
			}
			else cout << "Nema ga u stablu\n";
		}
		if (choice == 4) {
			string tmp;
			cout << " Unesi string: ";
			cin >> tmp;
			deleteNode(root, tmp);
		}if (choice == 5) {
			string tmp;
			cout << " Unesi string: ";
			cin >> tmp;
			cout << "Broj manjih od njega: ";
			cout << smallerKeys(root, tmp) << endl;
		}if (choice == 6) {
			cout << root;
			cout << endl;
		}if (choice == 7) {
			exit(1);
		}

	}
}