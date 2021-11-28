#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>

//TO DO PRELAMANJE U OPSTEM SLUCAJU, BRISANJE, SEARCH
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

	Node(const Node& n);
	Node(Node&& n);
	Node& operator=(const Node& n);
	Node& operator=(Node&& n);
	~Node();
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

//a b 

void deleteKeys(Node* nod) {
	for (int i = 0; i < nod->currElems; i++) {
		delete nod->data[i];
		nod->data[i] = nullptr;
	}
	//nod->currElems = 0;
}


void nodeSeparating(Node* separate, string d, int pos, int rl) {//pos nam pozicija prelomnog cvora, rl = 1 kad ima desnog brata
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
			Node* rootPointers[order + 2];
			for (int i = 0; i <= fath->currElems; i++)rootPointers[i] = fath->pointers[i], fath->pointers[i] = nullptr;
			rootPointers[fath->currElems + 1] = thirdNode;
			fath->currElems--;
			//pravimo 2 nova cvora, prelamamo koren
			insertNode(fath, nodeData[0]);
			insertNode(fath, nodeData[1]);
			Node* left = fath->pointers[0], *right = fath->pointers[1];
			for (int i = 0; i <= left->currElems; i++) {
				left->pointers[i] = rootPointers[i];
				left->pointers[i]->father = left;
			}
			for (int i = 0; i <= right->currElems; i++) {
				right->pointers[i] = rootPointers[i+left->currElems+1];
				right->pointers[i]->father = left;
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

Node* rightBrother(Node* tmp, int& pointerIndex, int& keyIndex) {
	Node* fath = tmp->father;
	for (int i = 0; i <= fath->currElems; i++) {
		if (fath->pointers[i] == tmp)pointerIndex = i;
	}
	keyIndex = pointerIndex - 1;
	if(pointerIndex +1 <= fath->currElems)return fath->pointers[pointerIndex + 1];
	return nullptr;
}

Node* leftBrother(Node* tmp, int& pointerIndex, int& keyIndex) {
	Node* fath = tmp->father;
	for (int i = 0; i <= fath->currElems; i++) {
		if (fath->pointers[i] == tmp)pointerIndex = i;
	}
	keyIndex = pointerIndex - 1;
	if(pointerIndex-1 >= 0)return fath->pointers[pointerIndex - 1];
	return nullptr;
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
				nodeSeparating(place, d, positionOfSon,0);
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
					nodeSeparating(place, d, pointerIndex, separatingrl);
				}
			}

		}
	}
}


int main() {
	Node* root = new Node;
	root->leaf = 0;
	root->root = 1;
	insertNode(root, "a");
	insertNode(root, "c");
	insertNode(root, "b");
	insertNode(root, "d");
	insertNode(root, "e");
	insertNode(root, "i");
	insertNode(root, "j");
	insertNode(root, "k");
	/*insertNode(root, "z");
	insertNode(root, "x");
	insertNode(root, "w");
	insertNode(root, "y");
	insertNode(root, "v");
	insertNode(root, "l");
	insertNode(root, "f");
	insertNode(root, "g");
	insertNode(root, "h");
	insertNode(root, "m");
	insertNode(root, "u");
	insertNode(root, "n");
	insertNode(root, "o");*/
	cout << root;
}