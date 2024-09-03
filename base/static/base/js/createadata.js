let filename_showplace = document.getElementById("showplace");
filename_showplace.classList.add("max-h-60", "overflow-y-auto");

// proccess button click
async function processbtn(event) {
  event.preventDefault();
  toggleLoading(true, "processbutton");
  const file_upload = document.getElementById("file_upload");
  const files = file_upload.files;

  if (files.length === 0) {
    alert("Please Upload Unless one File");
    toggleLoading(false, "processbutton");
    return;
  }

  const formData = new FormData();

  for (let i = 0; i < files.length; i++) {
    formData.append("files", files[i]);
  }
  const csrftoken = getCookie("csrftoken");
  const result = await fetchAPI("/api/upload", formData, csrftoken);

  toggleLoading(false, "processbutton");

  let adata_results = document.getElementById("adata_results");
  adata_results.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add(
    "max-w-md",
    "space-y-1",
    "text-gray-500",
    "list-disc",
    "list-inside"
  );
  result["adata_results"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  adata_results.appendChild(ul);
  let nextbtn = document.getElementById("nextbtn");
  nextbtn.classList.remove("hidden");
}

// file upload
document.getElementById("file_upload").addEventListener("change", function (e) {
  filename_showplace.innerHTML = "";
  files = e.target.files;
  if (files.length > 0) {
    for (let i = 0; i < files.length; i++) {
      const file = files[i];
      createFileDiv(file);
    }
  } else {
    filename_showplace.textContent = "No Data Upload Yet.";
  }
});

function createFileDiv(file) {
  const filediv = document.createElement("div");
  filediv.classList.add(
    "flex",
    "w-full",
    "justify-between",
    "items-center",
    "mb-2",
    "p-2",
    "bg-gray-100",
    "rounded"
  );

  const filename = document.createElement("div");
  filename.textContent = file.name;
  filename.classList.add("truncate", "w-1/2");

  const filesize = document.createElement("div");
  filesize.textContent = formatFileSize(file.size);
  filesize.classList.add("text-sm", "text-gray-500");

  const delfile = document.createElement("button");
  delfile.textContent = "X";
  delfile.classList.add(
    "bg-red-500",
    "text-white",
    "px-2",
    "py-1",
    "rounded",
    "hover:bg-red-600"
  );
  delfile.onclick = function () {
    filediv.remove();
    if (filename_showplace.children.length === 0) {
      filename_showplace.textContent = "No Data Upload Yet.";
    }
  };

  filediv.appendChild(filename);
  filediv.appendChild(filesize);
  filediv.appendChild(delfile);
  filename_showplace.appendChild(filediv);
}

function formatFileSize(bytes) {
  if (bytes < 1024) return bytes + " bytes";
  else if (bytes < 1048576) return (bytes / 1024).toFixed(2) + " KB";
  else if (bytes < 1073741824) return (bytes / 1048576).toFixed(2) + " MB";
  else return (bytes / 1073741824).toFixed(2) + " GB";
}
