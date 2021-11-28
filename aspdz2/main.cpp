#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>

/*allData.push_back(*separate->data[pos]); // prelomni kljuc iz oca
		//uzima desnog, ako nema desnog onda levog brata koji je takodje popunjen
		Node* Brother = (separate->father->pointers[pos + 1] ? separate->father->pointers[pos + 1] : separate->father->pointers[pos - 1]);
		for (int i = 0; i < Brother->currElems; i++)allData.push_back(*Brother->data[i]);
		Node* fath = separate->father;
		sort(allData.begin(), allData.end());

		deleteKeys(separate);
		deleteKeys(Brother);

		Node* firstNode = separate->father->pointers[pos + 1] ? separate : separate->father->pointers[pos - 1];

		for (int i = 0; i < firstNum; i++) {
			firstNode->data[i] = new string{ allData[i] };
		}
		inserted += firstNum + 1;

		Node* secondNode = separate->father->pointers[pos + 1] ? separate->father->pointers[pos + 1] : separate;
		for (int i = 0; i < secondNum; i++) {
			secondNode->data[i] = new string{ allData[i + firstNum + 1] };
		}
		inserted += secondNum;

		if (inserted != n) {// ako se uzme 3 granicnik mora da ima 3 cvora
			inserted++;

			Node* thirdNode = new Node;
			for (int i = 0; i < thirdNum; i++) {
				thirdNode->data[i] = new string{ allData[i + firstNum + secondNum + 2] };
			}
			inserted += thirdNum;
		}


		//napravimo 3 cvora i saljemo sad u oca
		if (separate->father->currElems + 2 < separate->father->maxPointers() - 1) {
			//ubacujemo ta dva kljuca i prepovezujemo
		}
		else {
			//cuvamo pokazivace iz oca pre pos, oca delimo
		}*/


//TO DO PRELAMANJE U OPSTEM SLUCAJU, BRISANJE, SEARCH
using namespace std;

const int nn = 5; //red stabla br pokazivaca
const int order = 5; //2 * floor((2 * nn - 2) / 3) + 1; OVO JE MAX BR STRINGOVA, SAMO KOREN IMA TOLIKO OSTALI
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

		while (!qu.empty()) {
			Node* tmp = qu.front();
			qu.pop();
			
			if (*tmp->data[0] == "XXXXX")os << endl;
			else {
				os << "[";
				for (int i = 0; i < tmp->currElems; i++) {
					os << *tmp->data[i] << " ";
				}
				os << "]   ";
				qu.push(gran);
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
	int greater = 0;
	while (temp) {
		find = temp;
		for (int i = 0; i < temp->currElems; i++) {
			if (*temp->data[i] > d) {
				temp = temp->pointers[i], pos = i;
				if (!temp)return find;
			}
			else {
				greater++;
			}
		}
		if (greater == temp->currElems)pos = temp->currElems, temp = temp->pointers[temp->currElems]; //ako je veci od svih
		greater = 0;
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
		int brotherPos = rl == 1 ? pos : pos - 1; //pos je pzicija pointera 
		Node* fath = separate->father;
		Node* brother = fath->pointers[brotherPos];
		for (int i = 0; i < separate->currElems; i++)allData.push_back(*separate->data[i]);
		for (int i = 0; i < brother->currElems; i++)allData.push_back(*brother->data[i]);
		allData.push_back(d);
		int keyInd = (pos == 0) ? 0 : pos - 1;
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
		delete fath->data[pos - 1];
		fath->data[keyInd] = nullptr;
		if (fath->currElems + nodeData.size() < fath->maxKeys()) {
			fath->data[keyInd] = new string{ nodeData[0] };
			fath->data[keyInd+1] = new string{ nodeData[1] };
			fath->pointers[pos + 1] = thirdNode;
			fath->currElems++;
		}
		else {
			cout << "nema mesta";
		}

	}

}

void nodePouring(Node* pour, Node* brother, int pos, string d, int rl) {//pos lokacija prelomnog
	vector<string> allData;
	Node* fath = pour->father;
	for (int i = 0; i < pour->currElems; i++)allData.push_back(*pour->data[i]);
	for (int i = 0; i < brother->currElems; i++)allData.push_back(*brother->data[i]);
	int pozition = (pos == 0 ? 0 : pos - 1);
	string s = *pour->father->data[pozition];
	allData.push_back(*pour->father->data[pozition]);
	allData.push_back(d);

	sort(allData.begin(), allData.end());

	int mid = allData.size() / 2;
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
				
				if (rightBrother(place, pointerIndex, keyIndex)) {//da li je u opsegu
					separatingrl = 1;
					if (!rightBrother(place, pointerIndex, keyIndex)->isFull()) {
						help = rightBrother(place, pointerIndex, keyIndex);
						flag = 1;
						rl = 1; 
					}
				}
				else if (leftBrother(place, pointerIndex, keyIndex)) { //da li je u opsegu
					if (!leftBrother(place, pointerIndex, keyIndex)->isFull()) {
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
	insertNode(root, "z");
	insertNode(root, "x");
	cout << root;
}