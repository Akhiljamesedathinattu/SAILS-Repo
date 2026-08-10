import torch
from torch import nn, device
from torch.utils.data import Dataset, DataLoader, Subset
from torchvision import datasets
from torchvision.transforms import v2


# Define neural network
class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()

        self.linear1 = nn.Linear(in_features=28 * 28, out_features=1024)
        self.relu1 = nn.ReLU
        self.linear2 = nn.Linear(in_features=1024, out_features=512)
        self.relu2 = nn.ReLU
        self.linear3 = nn.Linear(in_features=512, out_features=256)
        self.relu3 = nn.ReLU
        self.linear4 = nn.Linear(in_features=256, out_features=10)

    def forward(self, x):
        x = self.flatten(x)
        x = self.linear1(x)
        x = self.relu1(x)
        x = self.linear2(x)
        x = self.relu2(x)
        x = self.linear3(x)
        x = self.relu3(x)
        logits = self.linear4(x)
        out = logits
        return out


def load_data():
    # Download train dataset
    full_training_data = datasets.FashionMNIST(root='./data', train=True,download=True,transform=v2.Compose([v2.ToImage(), v2.ToDtype(torch.float32,scale=True)]))
    print(len(full_training_data))

    # Download test dataset
    full_test_data = datasets.FashionMNIST(root='./data', train=False,download=True,transform=v2.Compose([v2.ToImage(), v2.ToDtype(torch.float32,scale=True)]))
    print(len(full_test_data))

    return full_training_data, full_test_data


def train(mydataloader, model, loss_fn, optimizer, device, epochs):
    size = len(mydataloader.dataset)
    for epoch in range(epochs):
        model.train()
        for batch, (X, y) in enumerate(mydataloader):
            X, y = X.to(device), y.to(device)

            # Prediction error
            pred = model(X)
            loss = loss_fn(pred, y)

            # Backpropagation
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if batch % 100 == 0:
                loss, current = loss.item(), (batch + 1) * len(X)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")


def test(mydataloader, model, loss_fn, device):
    size = len(mydataloader.dataset)
    num_batches = len(mydataloader)
    model.eval()

    test_loss, correct = 0, 0
    with torch.no_grad():
        for X, y in mydataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) ==
                        y).type(torch.float).sum().item()

    test_loss /= num_batches
    correct /= size

    print(f"Test Error: \n Accuracy: {(100 * correct):>0.1f}%, Avg loss: {test_loss: > 8f} \n")


def main():
    # Load dataset
    training_data, test_data = load_data()

    batch_size = 64
    train_dataloader = DataLoader(training_data, batch_size=batch_size,
                                  shuffle=True)
    test_dataloader = DataLoader(test_data, batch_size=batch_size,
                                 shuffle=True)

    for X, y in train_dataloader:
        print(f"Shape of X [N, C, H, W]: {X.shape}")
        print(f"Shape of y: {y.shape} {y.dtype}")
        break

    # Load accelerators
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using {device} device")

    # Initialize the network
    model = NeuralNetwork().to(device)
    print(model)

    # Optimize
    loss_fn = nn.CrossEntropyLoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=1e-3)

    # Model training
    for epoch in range(10):
        train(train_dataloader, model, loss_fn, optimizer, device=device,
              epochs=1)
        test(test_dataloader, model, loss_fn, device)

    print('End')


if __name__ == "__main__":
    main()
